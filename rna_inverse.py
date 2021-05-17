import numpy as np
import pandas as pd
from nupack import *
import random
import statistics

def structure_to_bp(structure):
    """
    Converts dot parens structure to dictionary of base pairs
    
    Args:
        structure (str): dot-parens structure notation
    
    Returns:
        bp (dict): dictionary of base pair indices
    """
    open_parens = [] # List of open parentheses
    bp = {} # Dict of base pairs
    for i, x in enumerate(structure):
        if x == '(':
            open_parens.append(i)
        elif x == ')':
            open_index = open_parens.pop()
            bp[open_index] = i
            bp[i] = open_index
    return bp

def check_sequence_constraints(structure, constraints):
    """
    Checks the primary sequence constraints given by the user,
    converts user constraints to full list of constraints
    
    Args:
        structure (str): desired dot-parens structure notation
        constraints (dict): primary constraints {pos: base} form
            e.g. constraints[12] = 'A'

    Returns:
        new_constraints (dict): primary constraints converted
            to account for base pairing
    """
    new_constraints = constraints.copy()
    bases = ['A', 'G', 'U', 'C']
    pairs = {'A':'U', 'G':'C', 'U':'A', 'C':'G'}
    bp = structure_to_bp(structure)
    for position in constraints.keys():
        if position in bp:
            bp_position = bp[position]
            old_base = constraints[position]
            new_constraints[bp_position] = pairs[old_base]
    return new_constraints

def make_random_sequence(structure, constraints=None):
    """
    Takes a given structure and generates a random inverse 
    nucleotide sequence

    Args:
        structure (str): desired dot-parens structure notation
        constraints (dict): full list of primary constraints

    Returns:
        sequence (str): generated random RNA sequence
    """
    bases = ['A', 'G', 'U', 'C']
    pairs = {'A':'U', 'G':'C', 'U':'A', 'C':'G'}
    sequence = ['_'] * len(structure)
    bp = structure_to_bp(structure)
    for i, x in enumerate(structure):
        if constraints:
            if i in constraints.keys():
                sequence[i] = constraints[i]
        if x == '.':
            sequence[i] = random.choice(bases)
        elif x == '(':
            sequence[i] = random.choice(bases)
            close_index = bp[i]
            sequence[close_index] = pairs[sequence[i]] 
    sequence = ''.join(sequence)
    return sequence

def initialize_sequences(structure, number_of_sequences, constraints=None):
    """
    Generates sequences, evaluates based on mfe score and
    similarity score, ranks sequences.

    Args:
        structure (str): target structure
        number_of_sequences (int): number of sequences to generate and initialize
        constraints (dict): primary sequence constraints

    Returns:
        dataframe (DataFrame): sorted DataFrame of random sequences
        initialized_sequence (str): best ranked random sequence
    """
    random_sequences = []
    nupack_structures = []
    mfe_list = []
    mismatch_list = []
    similarity_scores = []
    for i in range(number_of_sequences):
        random_sequence = make_random_sequence(structure, constraints)
        random_sequences.append(random_sequence)
        nupack_structure, mfe = nupack_analyze_sequence(random_sequence)
        nupack_structures.append(nupack_structure)
        mfe_list.append(mfe)
        mismatch = structure_differences(structure, nupack_structure)
        mismatch_list.append(mismatch)
        similarity_scores.append((len(structure) - len(mismatch))/len(structure))
    dataframe = pd.DataFrame(list(zip(random_sequences, nupack_structures, mfe_list, mismatch_list, similarity_scores)), 
                             columns=['random_sequences', 'nupack_structures', 'mfe', 'mismatch', 'similarity_scores'])
    dataframe = dataframe.sort_values(['similarity_scores', 'mfe'], ascending=(False, True))
    initialized_sequence = dataframe['random_sequences'].iloc[0]
    return dataframe, initialized_sequence

def nupack_analyze_sequence(sequence):
    """
    Takes given sequence and runs it through Nupack package to
    return the lowest energy RNA secondary structure
    
    Args:
        sequence (str): RNA sequence
    
    Returns:
        nupack_structure (str): predicted RNA structure from Nupack
            package
        mfe (float): minimum free energy of sequence folding into structure
    """
    model = Model(material='rna', celsius=37) # Physical model
    A = Strand(sequence, name='A')
    tube = Tube({A: 1e-8}, complexes = SetSpec(max_size=1),
               name='Tube 1')
    results = tube_analysis(tubes=[tube], compute=['pairs', 'mfe'],
                            model=model)
    tube_result = results['(A)']
    nupack_structure = str(tube_result.mfe[0].structure)
    mfe = tube_result.mfe[0].energy
    return nupack_structure, mfe

def mutate_random_position(mismatch, sequence, structure, constraints=None):
    """
    Mutates a random position in sequence from list of indices
    in mismatch list
    
    Args:
        mismatch (list): list of indices from which random position
            can be mutated
        sequence (str): RNA sequence
        structure (str): dot-parens notation of desired RNA structure

    Returns:
        mutated_sequence (str): mutated RNA sequence
    """
    if constraints:
        for position in constraints.keys():
            mismatch = mismatch.remove(position) # Adding primary constraints
    bases = ['A', 'G', 'U', 'C']
    pairs = {'A':'U', 'G':'C', 'U':'A', 'C':'G'}
    bp = structure_to_bp(structure)
    mutate_position = random.choice(mismatch) # Question!
    sequence = list(sequence)
    sequence[mutate_position] = random.choice(bases)
    if mutate_position in bp:
        mutate_bp_position = bp[mutate_position]
        sequence[mutate_bp_position] = pairs[
            sequence[mutate_position]]
    mutated_sequence = ''.join(sequence)
    return mutated_sequence
        
def structure_differences(structure, nupack_structure):
    """
    Returns a list of indices that correspond to mismatched bases
    in ideal and predicted RNA structures

    Args:
        structure (str): target secondary structure
        nupack_structure (str): predicted structure from sequence

    Returns:
        mismatch (list): list of mismatched bases indices
    """
    mismatch = [] # List of mismatched bases indices
    structure = list(structure)
    nupack_structure = list(nupack_structure)
    for i in range(len(structure)):
        if structure[i] != nupack_structure[i]:
            mismatch.append(i)
    return mismatch

def compare_mutate_sequence(random_sequence, structure, nupack_structure, constraints=None):
    """
    Compares the ideal structure to the predicted structure of
    generated random RNA sequence from Nupack
    Mutates RNA sequence randomly if structures don't match
    
    Args:
        random_sequence (str): random RNA sequence from 
            random_sequence function
        structure (str): desired dot-parens structure notation
        nupack_structure (str): dot-parens structure of random sequence
            predicted from Nupack package

    Returns:
        mutated_mismatch (list): list of mismatched bases indices for mutated structure
        mutated_sequence (str): mutated sequence
    """
    mutated_mismatch = []
    mutated_sequence = random_sequence
    if structure != nupack_structure:
        mismatch = structure_differences(structure, nupack_structure)
        mutated_sequence = mutate_random_position(mismatch, random_sequence, structure)
        mutated_nupack_structure, _ = nupack_analyze_sequence(mutated_sequence)
        mutated_mismatch = structure_differences(structure, mutated_nupack_structure)
    return mutated_mismatch, mutated_sequence

def mutate_sequence_iterate(random_sequence, structure, nupack_structure, iterations, constraints=None):
    """
    Iteratively compares predicted structure from Nupack to ideal structure,
    mutates random sequence if mismatched, and adopts mutated sequence under
    certain conditions
    
    Args:
        random_sequence (str): starting sequence
        structure (str): target secondary structure
        nupack_structure (str): predicted secondary structure
        iterations (int): number of iterations/walks to perform

    Returns:
        result (dict): success, iteration to correct structure, final sequence, all sequences in order
    """
    mismatch = structure_differences(structure, nupack_structure)
    sequences = []
    for i in range(iterations):
        mutated_mismatch, mutated_sequence = compare_mutate_sequence(
            random_sequence, structure, nupack_structure)
        if len(mutated_mismatch) == 0: # Set stopping condition if structure matches
            random_sequence = mutated_sequence
            mismatch = mutated_mismatch
            break
        if len(mutated_mismatch) < len(mismatch):
            random_sequence = mutated_sequence
            mismatch = mutated_mismatch
            nupack_structure, _ = nupack_analyze_sequence(random_sequence)
        elif len(mutated_mismatch) == len(mismatch):
            random_sequence = mutated_sequence
            mismatch = mutated_mismatch
            nupack_structure, _ = nupack_analyze_sequence(random_sequence)
        elif len(mutated_mismatch) > len(mismatch):
            random_number = random.random()
            if random_number >= 0.95:
                random_sequence = mutated_sequence
                mismatch = mutated_mismatch
                nupack_structure, _ = nupack_analyze_sequence(random_sequence)
        sequences.append(random_sequence)
    result = {}    
    if len(mismatch) > 0:
        result['success'] = False
    else:
        result['success'] = True
    result['iteration'] = i
    result['final_sequence'] = random_sequence
    result['all_sequences'] = sequences
    return result