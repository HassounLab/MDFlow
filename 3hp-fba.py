import cobra.test
from cobra import Model, Reaction, Metabolite
from cobra.io import load_json_model
from cobra.flux_analysis.moma import add_moma
import json
import re
import numpy as np
import random
import sys
from hashlib import md5

if len(sys.argv) != 5:
    print('Usage: %s <pathway: p1, p2> <s1, s2, fixed> <fba, moma> <seed>' % sys.argv[0])
    exit()
assert(sys.argv[1] in ['p1', 'p2'])
isPathway1 = sys.argv[1] == 'p1'
mode = sys.argv[2]
method = sys.argv[3]
seed = int(sys.argv[4])
assert(mode in ['s1', 's2', 'fixed'])
assert(method in ['fba', 'moma'])

random.seed(seed)
np.random.seed(seed)

COUPLING_MEAN = 0.01
EXPERIMENTS = 10000

print('Command:', sys.argv)

print('Coupling: %f mean' % COUPLING_MEAN)
print('Total experiments: %d' % EXPERIMENTS)

# Import reactions and filter out duplicates

def reactionKey(rxn):
    return '='.join(sorted(['+'.join(c for c, f in rxn[0]), # direction-agnostic
                            '+'.join(c for c, f in rxn[1])]))
def reactionRowKey(key):
    return md5(key.encode('ascii')).hexdigest()[0:4]

blockedKeys = set([
    '3318', '2571', 'aabc', 'e541', '6dbd', '097b',
    '374a', '4774', '7b33', '0e18', '5d0a', '4b66',
    '8f1e', 'c9a4', 'f234', '251c', '8056', '84c7',
    '6c29', 'dce1', '6281', '5f86', 'e3aa', 'dfd7',
    '5c8f'
])

def importReactions(pathway, scenario):
    reactionEnzymes = {} # key: ec list

    allReactions = []
    with open('data/reactions-3hp.tsv', 'r') as f:
        for line in f:
            label, substrates, enzyme, reactions, products, equationsJson = line.strip().split('\t')
            if label != ('3hp-p%d-s%d' % (pathway, scenario)) and \
                    label != ('p%d-s%d' % (pathway, scenario)):
                continue # skip irrelevant pathway or scenario
            newReactions = [rxn for rxn in json.loads(equationsJson) if reactionRowKey(reactionKey(rxn)) not in blockedKeys]
            for rxn in newReactions:
                k = reactionKey(rxn)
                if k not in reactionEnzymes:
                    reactionEnzymes[k] = []
                reactionEnzymes[k].append(enzyme)
            allReactions.extend(newReactions)

    reactionsByKey = {}
    for rxn in allReactions:
        reactionsByKey[reactionKey(rxn)] = rxn
    reactionKeys = sorted(list(reactionsByKey.keys()))
    reactions = [reactionsByKey[key] for key in reactionKeys]

    reactionOrder = [0] * len(allReactions) # index in the json: index in reactions array
    for i, rxn in enumerate(allReactions):
        reactionOrder[i] = reactionKeys.index(reactionKey(rxn))

    print('%s: Imported %d reactions, %d unique found' %
          (scenario, len(allReactions), len(reactions)))

    return reactions, reactionOrder, reactionEnzymes

reactionsS1, reactionOrderS1, reactionEnzymesS1 = importReactions(1 if isPathway1 else 2, 1)
reactionsS2, reactionOrderS2, reactionEnzymesS2 = importReactions(1 if isPathway1 else 2, 2)

model = load_json_model('data/iML1515.json')

# Setup M9 growth medium

def countAtoms(formula):
    atoms = {}
    for atom, count in re.findall(r'([A-Z]{1}[a-z]*)([0-9]*)', formula):
        if atom not in atoms:
            atoms[atom] = 0
        if count == '':
            atoms[atom] += 1
        else:
            atoms[atom] += int(count)
    return atoms
assert(countAtoms('C6H13NO2') == {'C': 6, 'H': 13, 'N': 1, 'O': 2})
assert(countAtoms('C6H13PbO2') == {'C': 6, 'H': 13, 'Pb': 1, 'O': 2})
assert(countAtoms('C6H13Pb51O2') == {'C': 6, 'H': 13, 'Pb': 51, 'O': 2})

for rxn in model.reactions:
    if not rxn.id.startswith('EX_'): # not exchange reaction
        continue

    assert(len(rxn.metabolites) == 1)
    formula = next(iter(rxn.metabolites.keys())).formula
    atoms = countAtoms(formula)

    if 'C' in atoms:
        rxn.lower_bound = 0 # suppress carbon sources

model.reactions.get_by_id('EX_glc__D_e').lower_bound= -10
model.reactions.get_by_id('EX_o2_e').lower_bound = -20

# Set lower bound for biomass production to 10% of normal

biomassCoreRxn = model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M')
biomassWtRxn = model.reactions.get_by_id('BIOMASS_Ec_iML1515_WT_75p37M')

biomassCoreRxn.objective_coefficient = 1.0
biomassWtRxn.objective_coefficient = 0
model.optimize()
biomassCoreFlux = biomassCoreRxn.flux
print('Normal core biomass growth', biomassCoreFlux)

biomassCoreRxn.objective_coefficient = 0
biomassWtRxn.objective_coefficient = 1.0
model.optimize()
biomassWtFlux = biomassWtRxn.flux
print('Normal wild type biomass growth', biomassWtFlux)

biomassCoreRxn.objective_coefficient = 0
biomassWtRxn.objective_coefficient = 0
biomassCoreRxn.lower_bound = biomassCoreFlux * 0.1
biomassWtRxn.lower_bound = biomassWtFlux * 0.1

# Create 3-HP sink reaction as a target for optimization and determine fluxes under native conditions

targetRxn = Reaction('3HP_SINK')
targetRxn.name = 'Synthetic 3-HP Sink'
targetRxn.lower_bound = 0
targetRxn.upper_bound = 1000
targetRxn.add_metabolites({model.metabolites.get_by_id('3hpp_c'): -1})
model.add_reactions([targetRxn])

targetRxn.objective_coefficient = 1

# Add 3-HP synthesis steps (any new steps will be disabled initially)

pathwayRxns = []
if isPathway1:
    rxn = Reaction('3HP_P1_1')
    rxn.name = 'Synthetic 3-HP Pathway 1 Step 1'
    rxn.add_metabolites({ # Balanced
        model.metabolites.get_by_id('malcoa_c'): -1, # C24H33N7O19P3S
        model.metabolites.get_by_id('h_c'): -1, # H
        model.metabolites.get_by_id('nadph_c'): -1, # C21H26N7O17P3
        model.metabolites.get_by_id('nadp_c'): 1, # C21H25N7O17P3
        model.metabolites.get_by_id('coa_c'): 1, # C21H32N7O16P3S
        model.metabolites.get_by_id('msa_c'): 1, # C3H3O3
    })
    model.add_reactions([rxn])
    pathwayRxns.append(rxn)

    rxn.lower_bound = 0
    rxn.upper_bound = 0
    rxn.objective_coefficient = 0

    # 2nd step: already exists as MSAR
    pathwayRxns.append(model.reactions.get_by_id('MSAR'))
else:
    # 1st step: already exists as ASP1DC
    pathwayRxns.append(model.reactions.get_by_id('ASP1DC'))

    rxn = Reaction('3HP_P2_2')
    rxn.name = 'Synthetic 3-HP Pathway 2 Step 2'
    rxn.add_metabolites({
        model.metabolites.get_by_id('ala_B_c'): -1, # beta-alanine, C3H7NO2
        model.metabolites.get_by_id('akg_c'): -1, # 2-oxoglutarate, C5H4O5
        model.metabolites.get_by_id('glu__L_c'): 1, # glutamate, C5H8NO4
        model.metabolites.get_by_id('msa_c'): 1, # malonate semialdehyde, C3H3O3
    })
    model.add_reactions([rxn])
    pathwayRxns.append(rxn)

    rxn.lower_bound = 0
    rxn.upper_bound = 0
    rxn.objective_coefficient = 0

    # 3rd step: already exists as MSAR
    pathwayRxns.append(model.reactions.get_by_id('MSAR'))

# Add developed reactions and sinks for unknown metabolites (all disabled)

def addDevelopedReactions(model, scenario, reactions):
    developedRxns = {}

    for i, (substrates, products) in enumerate(reactions):
        if len(substrates) == 0:
            assert(len(products) == 0)
            continue

        rxn = Reaction('NEW_%s_%d' % (scenario, i))
        rxn.name = 'Developed %s reaction #%d' % (scenario, i)

        for metId in set(c for c, f in substrates) | set(c for c, f in products):
            try:
                met = model.metabolites.get_by_id(metId)
            except KeyError:
                met = Metabolite(metId, compartment='c')
                model.add_metabolites([met])

                # Any new produced metabolites require a sink, so add it in

                sinkId = metId + '_SINK'
                try:
                    sink = model.reactions.get_by_id(sinkId)
                except KeyError:
                    sink = Reaction(sinkId)
                    sink.name = metId + ' Sink'
                    sink.lower_bound = 0 # not a source!
                    sink.upper_bound = 1000
                    sink.add_metabolites({met: -1})
                    model.add_reactions([sink])

        metabolites = {}
        for cpd, factor in substrates:
            metabolites[model.metabolites.get_by_id(cpd)] = -factor
        for cpd, factor in products:
            metabolites[model.metabolites.get_by_id(cpd)] = factor
        rxn.add_metabolites(metabolites)

        model.add_reactions([rxn])
        developedRxns[i] = rxn

        rxn.lower_bound = 0
        rxn.upper_bound = 0
        rxn.objective_coefficient = 0

    return developedRxns

developedRxnsS1 = addDevelopedReactions(model, 's1', reactionsS1)
developedRxnsS2 = addDevelopedReactions(model, 's2', reactionsS2)

# Compute baseline behavior

print('3-HP Production, Baseline')
baseline = model.optimize()
baselineTargetFlux = baseline['3HP_SINK']
print('    Growth %f (core), %f (wt)' % (baseline[biomassCoreRxn.id], baseline[biomassWtRxn.id]))
print('    Glucose %f -> 3HP %f' % (baseline['EX_glc__D_e'], baselineTargetFlux))

# Set up MOMA constraints if needed

if method == 'moma':
    add_moma(model, solution=baseline)

# Compute engineered behavior with pathway reactions enabled

for rxn in pathwayRxns:
    rxn.lower_bound = -1000
    rxn.upper_bound = 1000

print('3-HP Production, w/Synthetic Path')
engineered = model.optimize()
normalTargetFlux = engineered['3HP_SINK']
print('    Growth %f (core), %f (wt)' % (engineered[biomassCoreRxn.id], engineered[biomassWtRxn.id]))
print('    Glucose %f -> 3HP %f' % (engineered['EX_glc__D_e'], normalTargetFlux))

if abs(normalTargetFlux - baselineTargetFlux) < 0.01:
    print('The addition of the synthesis pathway resulted in no substantial increase of yield!')
    print('Make sure the model and its configuration is appropriate for this experiment.')
    exit()

# Compute reaction flux bounds that maximize disruption (assuming certain coupling %)

def computeFluxBounds(direction, minFlux, maxFlux, coupling):
    if direction < 0: # reverse
        return (-1000, maxFlux + coupling * (minFlux - maxFlux))
    elif direction > 0: # forward
        return (minFlux + coupling * (maxFlux - minFlux), 1000)
    else: # bidirectional
        return  (-1000, 1000)

def determineFluxLimits(model, developedRxns, reactionOrder):
    rxnFluxLimits = {} # rxn index: (direction, min flux, max flux)
    for i, rxn in developedRxns.items():
        with model as model:
            targetRxn.objective_coefficient = 0

            # Find min forward flux
            rxn.objective_coefficient = -1
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            solution = model.optimize()
            minFwdFlux = solution[rxn.id]

            # Find max forward flux
            rxn.objective_coefficient = 1
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            solution = model.optimize()
            maxFwdFlux = solution[rxn.id]

            # Find min reverse flux
            rxn.objective_coefficient = -1
            rxn.lower_bound = -1000
            rxn.upper_bound = 0
            solution = model.optimize()
            minRevFlux = solution[rxn.id]

            # Find max reverse flux
            rxn.objective_coefficient = 1
            rxn.lower_bound = -1000
            rxn.upper_bound = 0
            solution = model.optimize()
            maxRevFlux = solution[rxn.id]

            targetRxn.objective_coefficient = 1
            rxn.objective_coefficient = 0

            # Compute coupled flux and disruption in forward direction
            rxn.lower_bound, rxn.upper_bound = \
                computeFluxBounds(+1, minFwdFlux, maxFwdFlux, COUPLING_MEAN)
            solution = model.optimize()
            coupledFwdFlux = solution[rxn.id]
            targetFwdFlux = solution[targetRxn.id]
            disruptionFwd = (targetFwdFlux - normalTargetFlux) / normalTargetFlux

            # Compute coupled flux and disruption in reverse direction
            rxn.lower_bound, rxn.upper_bound = \
                computeFluxBounds(-1, minRevFlux, maxRevFlux, COUPLING_MEAN)
            solution = model.optimize()
            coupledRevFlux = solution[rxn.id]
            targetRevFlux = solution[targetRxn.id]
            disruptionRev = (targetRevFlux - normalTargetFlux) / normalTargetFlux

            # record reaction bounds from the most disruptive case
            # (negative % preferred, and the more negative the better)
            if disruptionFwd < disruptionRev:
                print('    Rxn %3d ( ->, %5d %5d): %f disruption\t%s' %
                      (i, minFwdFlux, maxFwdFlux, disruptionFwd, str(rxn)))
                rxnFluxLimits[i] = (+1, minFwdFlux, maxFwdFlux)
            elif disruptionRev < disruptionFwd:
                print('    Rxn %3d (<- , %5d %5d): %f disruption\t%s' %
                      (i, minRevFlux, maxRevFlux, disruptionRev, str(rxn)))
                rxnFluxLimits[i] = (-1, minRevFlux, maxRevFlux)
            else:
                print('    Rxn %3d (<->, %5d %5d): %f disruption\t%s' %
                      (i, -1000, 1000, disruptionFwd, str(rxn)))
                rxnFluxLimits[i] = (0, -1000, 1000)
    return rxnFluxLimits

print('Developed S1 reaction flux bounds:')
rxnFluxLimitsS1 = determineFluxLimits(model, developedRxnsS1, reactionOrderS1)
print('Developed S2 reaction flux bounds:')
rxnFluxLimitsS2 = determineFluxLimits(model, developedRxnsS2, reactionOrderS2)

# Perform coupling percentage sampling experiments

def runCouplingVsDisruptionSimulations(model, filename,
    developedRxnsA, rxnFluxLimitsA, activeRxnNumA,
    developedRxnsB=None, rxnFluxLimitsB=None, activeRxnNumB=None):
    couplings = []
    disruptions = []

    for i in range(EXPERIMENTS):
        if i % 100 == 0:
            print('%s: experiment #%d...' % (filename, i))

        couplingMean = None
        while couplingMean is None or couplingMean <= 0:
            if couplingMean is not None:
                print('    Resampling coupling mean: previous was %f' % couplingMean)
            couplingMean = np.random.normal(COUPLING_MEAN, COUPLING_MEAN / 3)

        activeRxnIndicesA = random.sample(developedRxnsA.keys(), activeRxnNumA)
        assert(len(activeRxnIndicesA) > 0)

        if developedRxnsB is not None:
            activeRxnIndicesB = random.sample(developedRxnsB.keys(), activeRxnNumB)
        else:
            activeRxnIndicesB = []

        with model as model:
            for j in activeRxnIndicesA:
                rxn = developedRxnsA[j]
                direction, minFlux, maxFlux = rxnFluxLimitsA[j]
                coupling = np.random.normal(couplingMean, couplingMean / 3) # 99.7% mass above 0
                rxn.lower_bound, rxn.upper_bound = \
                    computeFluxBounds(direction, minFlux, maxFlux, coupling)

            for j in activeRxnIndicesB:
                rxn = developedRxnsB[j]
                direction, minFlux, maxFlux = rxnFluxLimitsB[j]
                coupling = np.random.normal(couplingMean, couplingMean / 3)
                rxn.lower_bound, rxn.upper_bound = \
                    computeFluxBounds(direction, minFlux, maxFlux, coupling)

            targetRxn.objective_coefficient = 1
            disrupted = model.optimize()

            disruption = (disrupted[targetRxn.id] - normalTargetFlux) / normalTargetFlux

        couplings.append(couplingMean)
        disruptions.append(disruption)

    with open('output/' + filename, 'w') as f:
        json.dump({'coupling': couplings, 'disruption': disruptions}, f)

ACTIVE_RXN_FRACTIONS = [0.10, 0.25, 0.50]
ACTIVE_RXN_NUMBER = 5
if mode == 's1': # experiments selecting a % of active reactions in s1
    for activeRxnFraction in ACTIVE_RXN_FRACTIONS:
        runCouplingVsDisruptionSimulations(model,
                                           '3hp-p%d-seed%d-s1-%d%%.json' % (1 if isPathway1 else 2,
                                                    seed, activeRxnFraction * 100),
                                           developedRxnsS1, rxnFluxLimitsS1,
                                           round(len(developedRxnsS1) * activeRxnFraction))
elif mode == 's2': # experiments selecting a % of active reactions in s2
    for activeRxnFraction in ACTIVE_RXN_FRACTIONS:
        runCouplingVsDisruptionSimulations(model,
                                           '3hp-p%d-seed%d-s2-%d%%.json' % (1 if isPathway1 else 2,
                                                    seed, activeRxnFraction * 100),
                                           developedRxnsS2, rxnFluxLimitsS2,
                                           round(len(developedRxnsS2) * activeRxnFraction))
elif mode == 'fixed': # experiments selecting a fixed number of active reactions, from each and both
    runCouplingVsDisruptionSimulations(model,
                                       '3hp-p%d-seed%d-s1-%d.json' % (1 if isPathway1 else 2,
                                                seed, ACTIVE_RXN_NUMBER * 2),
                                       developedRxnsS1, rxnFluxLimitsS1, ACTIVE_RXN_NUMBER * 2)
    runCouplingVsDisruptionSimulations(model,
                                       '3hp-p%d-seed%d-s2-%d.json' % (1 if isPathway1 else 2,
                                                seed, ACTIVE_RXN_NUMBER * 2),
                                       developedRxnsS2, rxnFluxLimitsS2, ACTIVE_RXN_NUMBER * 2)
    runCouplingVsDisruptionSimulations(model,
                                       '3hp-p%d-seed%d-both-%d.json' % (1 if isPathway1 else 2,
                                                seed, ACTIVE_RXN_NUMBER),
                                       developedRxnsS1, rxnFluxLimitsS1, ACTIVE_RXN_NUMBER,
                                       developedRxnsS2, rxnFluxLimitsS2, ACTIVE_RXN_NUMBER)
