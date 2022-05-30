with open('wilmer.csv') as f:
    lines = f.read().split('\n')
    lines = [l.split('\t') for l in lines] 
    lines = [l for l in lines if len(l)==4]
    reactions = [l for l in lines if l[2].startswith('[c]')]
    reactions = [r[2].split(':')[1].strip() for r in reactions]
    forwards = [r.split(' --> ') for r in reactions if '-->' in r]
    twoways  = [r.split(' <==> ') for r in reactions if '<==>' in r]

def get_reagents_dict(reagents_string):
    reagents = reagents_string.split(' + ')
    reagents = [r.split() for r in reagents] 
    reagents = [['(1)', r[0]] if len(r)==1 else r for r in reagents] 

    reagents = [(reagent_name, float(multiplicity[1:-1]))
                for multiplicity,reagent_name in reagents]

    for (rnm, mlt) in reagents:
        if mlt != int(mlt): return None

    return {
            rnm:int(mlt) for rnm, mlt in reagents
    }
    return reagents_dict

forwards = [(get_reagents_dict(lhs),
             get_reagents_dict(rhs)) for lhs,rhs in forwards]
forwards = [(lhs,rhs) for lhs,rhs in forwards if lhs is not None and rhs is not None]

all_species = set(rnm for lhs,rhs in forwards
                      for rnm in list(lhs.keys())+list(rhs.keys()))
species_from_symbol = {':A'+str(i):rnm for i,rnm in enumerate(all_species)} 
symbol_from_species = {v:k for k,v in species_from_symbol.items()}


print('detected {} many species'.format(len(all_species)))

template = '''
GRAPH_NAME = LabelledPetriNet([SPECIES_SYMBOLS], 
    REACTION_DICTIONARY
) 
'''
def to_julia(list_of_reactions, graph_name = 'Brusselator'):
    relevant_species = {rnm for lhs,rhs in list_of_reactions
                        for rnm in list(lhs.keys())+list(rhs.keys())} 
    julia_code = template
    julia_code = julia_code.replace('GRAPH_NAME', 'Brusselator')
    julia_code = julia_code.replace('SPECIES_SYMBOLS', ', '.join(relevant_species))
    julia_code = julia_code.replace(
    'REACTION_DICTIONARY', ',\n    '.join(
       '{} => (({}) => ({}))'.format(':t'+str(i),
           ', '.join(symbol_from_species[rnm] for rnm, mult in lhs.items() for j in range(mult)),
           ', '.join(symbol_from_species[rnm] for rnm, mult in rhs.items() for j in range(mult))
       ) 
       for i,(lhs,rhs) in enumerate(list_of_reactions)))
    return julia_code

for i,reaction in enumerate(forwards):
    with open('junk/reaction-'+str(i)+'.jl','w') as f:
        f.write(to_julia([reaction]))

with open('reactions-combined.jl','w') as f:
    f.write(to_julia(forwards))

