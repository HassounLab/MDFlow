import cobra.test
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model, load_json_model
from cobra.flux_analysis import moma, add_moma
import json, sys

try:
    assert(len(sys.argv) == 2)
    method = sys.argv[1].lower()
    assert(method in ['fba', 'moma'])
except:
    print(('Usage: %s <fba, moma>') % sys.argv[0])
    exit()

model = load_json_model('data/iML1428-iso_Glucose.json')

# Set constrants for aerobic growth in glucose minimal media
model.reactions.get_by_id('EX_glc__D_e').lower_bound= -10
model.reactions.get_by_id('EX_o2_e').lower_bound = -15

# Load experiments and reaction data

experiments = [
    'Isozyme overexpression',
    ('glyA', 'ltaE'),
    ('ilvA', 'tdcB'),
    ('ilvE', 'avtA'),
    ('metC', 'malY'),

    'Substrate ambiguity',
    ('glnA', 'asnB'),
    ('pdxB', 'tdh'),
    ('serB', 'gph'),
    ('serB', 'hisB'),
    ('serB', 'ytjC'),

    'Ambiguity in metabolite transport',
    ('fes', 'setB'),
    ('ilvA', 'emrD'),
    ('ppc', 'ecfM'),
    ('ppc', 'yccT'),
    ('ptsI', 'fucP'),
    ('ptsI', 'xylE'),

    'Catalytic promiscuity',
    ('fes', 'thiL'),
    ('metC', 'alr'),
    ('pabA', 'menF'),
    ('pabB', 'menF'),
    ('pdxB', 'purF'),
    ('purK', 'anmK'),

    'Maintaining metabolic flux',
    ('ilvD', 'avtA'),
    ('purK', 'purE'),

    'Metabolic pathway bypass and/or alternate source of downstream substrate',
    ('carA', 'carB'),
    ('glyA', 'tdh'),
    ('glyA', 'yneH'),
    ('hisH', 'hisF'),
    ('pabA', 'pabB'),
    ('ptsI', 'galE'),
    ('serA', 'yneH'),
    ('serC', 'yneH'),

    'Regulatory effects',
    ('carA', 'ygiT'),
    ('glyA', 'rsd'),
    ('metC', 'fimE'),
    ('metR', 'metE'),
    ('purK', 'cpdA'),

    'Unknown',
    ('carA', 'cho'),
    ('carA', 'yncK'),
    ('purF', 'chbA'),
    ('yhhK', 'murI'),
    ('yhhK', 'purF'),
]

genes = { # from uniprot
    'yhhK': 'b3459',
    'yneH': 'b1524'
}
with open('data/kegg-list-eco.txt') as f:
    for line in f:
        b, name, comment = line.strip().split(maxsplit=2)
        genes[name.rstrip(';')] = b.split(':')[1]

enzymes = { # from uniprot
    'b3117': ['4.3.1.17']
}
with open('data/kegg-link-eco-ec.txt') as f:
    for line in f:
        ec, b = line.strip().split()
        b = b.split(':')[1]
        if b not in enzymes:
            enzymes[b] = []
        enzymes[b].append(ec.split(':')[1])

reactions = {}
with open('data/reactions-patrick.tsv', 'r') as f:
    for line in f:
        label, substrates, enzyme, reactionIds, products, equationsJson = line.strip().split('\t')
        assert(label.startswith('patrick'))
        if enzyme not in reactions:
            reactions[enzyme] = []
        reactions[enzyme].extend(json.loads(equationsJson))

# Determine baseline growth of the model

baseline = model.optimize()
normalGrowth = baseline.objective_value

normalReactionFluxes = {}
for rxn in model.reactions:
    normalReactionFluxes[rxn.id] = rxn.flux

# Set up MOMA constraints if needed

if method == 'moma':
    add_moma(model, solution=baseline)

# Go through all MS experiments

out = open('output/patrick-%s.tsv' % method, 'w')
out.write('Deletion\tGrowth Rate\tDeficient Reactions\tNormal Def. Reaction Flux\t' +
          'MS\tMS Enzymes\tCompensating Reaction\tComp. Rxn. Flux\tGrowth Rate with Comp. Rxn.\n\n')

for entry in experiments:
    if isinstance(entry, str):
        out.write('*** %s\n\n' % entry)
        continue
    else:
        knockout, replacer = entry

    with model as model:
        # Knock out the essential gene and measure growth

        koGene = genes[knockout]
        try:
            koReactions = cobra.manipulation.find_gene_knockout_reactions(model, [koGene])
        except KeyError:
            out.write(knockout + '\tNot in model\n\n')
            continue
        for rxn in koReactions:
            rxn.knock_out()

        koGrowth = model.slim_optimize()
        if method == 'moma':
            # MOMA replaces the objective function so we have to retrieve growth differently.
            # See cobra.flux_analysis.deletion._get_growth() and "moma_old_objective" in docs.
            koGrowth = model.solver.variables.moma_old_objective.primal

        # Retrieve replacer reactions and enzymes (sort everything for repeatability)

        deficientResults = [(rxn.build_reaction_string(True), '%.6f' % normalReactionFluxes[rxn.id])
                            for rxn in sorted(koReactions, key=lambda rxn: rxn.id)]

        try:
            replacerEnzymes = enzymes[genes[replacer]]
        except KeyError:
            replacerEnzymes = []

        replacerReactions = [eqn for ec, equations in reactions.items() if ec in replacerEnzymes
                                    for eqn in equations]
        replacerReactions.sort(key=lambda equation: str(equation))

        # Generate replacer reactions (all of them will be disabled initially)

        externalMetIds = set()
        modelReactionsKeys = set()
        modelReactions = []
        for i, (substrates, products) in enumerate(replacerReactions):
            rxn = Reaction('NEW_%d' % i)
            rxn.name = 'Developed reaction #%d' % i
            rxn.lower_bound = 0
            rxn.upper_bound = 0

            for metId in set(c for c, f in substrates) | set(c for c, f in products):
                assert(metId)

                try:
                    met = model.metabolites.get_by_id(metId)
                except KeyError:
                    met = Metabolite(metId, compartment='c')
                    met.name = 'External ' + metId
                    model.add_metabolites([met])
                    externalMetIds.add(metId)

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

            key = rxn.build_reaction_string()
            for m in externalMetIds:
                key = key.replace(m, 'EXTERNAL')
            assert(' --> ' in key)
            key = ' <=> '.join(sorted(key.split(' --> ')))

            if key not in modelReactionsKeys:
                modelReactions.append(rxn)
                modelReactionsKeys.add(key)
            else:
                continue
        model.add_reactions(modelReactions)

        # Enable each reaction individually and measure growth rate

        rescueResults = []
        rescuedGrowth = []
        for rxn in modelReactions:
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

            growth = model.slim_optimize()
            if method == 'moma':
                growth = model.solver.variables.moma_old_objective.primal

            rescueResults.append((rxn.build_reaction_string(True), '%.6f' % rxn.flux))
            rescuedGrowth.append(growth)

            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # Print out all information we've got for this experiment

        rowCount = max(len(deficientResults), len(replacerEnzymes), len(rescueResults))
        assert(len(rescueResults) == len(rescuedGrowth))
        for i in range(rowCount):
            if i == 0:
                out.write('\t'.join([
                        knockout,
                        '%.2f%%' % (100 * koGrowth / normalGrowth)]) + '\t')
            else:
                out.write('\t'.join([''] * 2) + '\t')

            if i < len(deficientResults):
                out.write('\t'.join(str(x) for x in deficientResults[i]) + '\t')
            else:
                out.write('\t'.join([''] * len(deficientResults[0])) + '\t')

            if i == 0:
                out.write(replacer + '\t')
            else:
                out.write('\t')

            if i < len(replacerEnzymes):
                out.write(replacerEnzymes[i] + '\t')
            else:
                out.write('\t')

            if i < len(rescueResults):
                out.write('\t'.join(str(x) for x in rescueResults[i]) + '\t')
                out.write('%.2f%%' % (100 * rescuedGrowth[i] / normalGrowth) + '\n')
            else:
                out.write('\n')

        out.write('\n')

out.close()
