import uuid
from collections import Counter

def generate_id():
    return str(uuid.uuid4())

def generate_alL_schemas(nbits: int) -> list:
    """Generates a list of all wildcard schemas of length `nbits`.

    Args:
        nbits (int): number of bits in the binary string

    Returns:
        list: all possible combinations of 0,1,* schemas of length `nbits`
    """
    alphabets = ["0", "1", "*"]
    # generate list of all schemas
    all_schemas = [''.join(i) for i in itertools.product(alphabets, repeat = nbits)]
    # subset of schemas containing wildcards (avoid duplication)
    wildcard_schemas = [s for s in all_schemas if '*' in s]
    # remove the schema containing only wildcards '***..'
    wildcard_schemas = wildcard_schemas[:-1]
    return wildcard_schemas

def generate_schemas_from_solution(bitstring, original=None, index=0) -> list:
    """Given a bitstring solution, this recursive function returns all wildcard allele schemas present in that solution.

    Args:
        solution (str): a bitstring solution, eg "010011011".
    
    Returns:
        list: a list of all wildcard allele schemas present in the solution.
    """
    if original is None:
        original = bitstring
    if index == len(bitstring):
        if bitstring == original:
            return []
        else:
            return [bitstring]
    else:
        combinations = []
        combinations += generate_schemas_from_solution(bitstring, original, index + 1)
        bitstring = bitstring[:index] + '*' + bitstring[index+1:]
        combinations += generate_schemas_from_solution(bitstring, original, index + 1)
        return combinations
    # schemata = []
    # for i, bit in enumerate(solution):
    #     # if this is the first bit, initialize the schemata list
    #     if i == 0:
    #         schemata.append(bit)
    #         schemata.append("*")
    #     else:
def generate_schemas_from_population(population:list) -> list:
    """Given a population of `nbit` length bitstring solutions, this returns all wildcard allele schemas present in that population.

    Args:
        population (list): a list of `nbit` length bitstring solutions.

    Returns:
        list: a list of all wildcard allele schemas present in the population.
    """
    if len(set([len(s) for s in population])) != 1:
        raise Exception("Population contains solutions of different lengths. All solutions must be the same length bitstring.")
    
    schemata = []
    for solution in population:
        schemata_matches = generate_schemas_from_solution(solution)
        schemata += schemata_matches
    return schemata

def count_schemata(schemata: list, counts_dict: dict) -> dict:
    for element in schemata:
        if element in counts_dict.keys():
            counts_dict[element] += 1
        else:
            counts_dict[element] = 1
    return counts_dict    

# def generate_schemas_from_population(population:list, schemata_set:list) -> list:
#     """Given a population of `nbit` length bitstring solutions, this returns all wildcard allele schemas present in that population.

#     Args:
#         population (list): a list of `nbit` length bitstring solutions.

#     Returns:
#         list: a list of all wildcard allele schemas present in the population.
#     """
#     if len(set([len(s) for s in population])) != 1:
#         raise Exception("Population contains solutions of different lengths. All solutions must be the same length bitstring.")
    
#     for solution in population:
#         schemata_matches = [s for s in schemata_set if fnmatch.fnmatch(solution, s)]
#         print("solution:", solution)
#         print("wildcard matches:", schemata_matches)
#         # wildcard_schemata = 0
#     pass
# def order(allele:str) -> int:
#     # number of non wildcard bits
    return sum(1 for bit in allele if bit != "*")

def length(schema:str) -> int:
    # find first non-wildcard
    for i, bit in enumerate(schema):
        if bit != "*":
            start = i
            break
    # find last non-wildcard
    for i, bit in enumerate(reversed(schema)):
        if bit != "*":
            end = (len(schema)-1)-i
            break
        
    # start = next(bit for bit in allele if bit != "*")
    # [*,0,0,*,1,*]
    # [^*]
    # end = next(bit for bit in reversed(allele) if bit != "*")
    l = end - start
    return l