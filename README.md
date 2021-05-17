# CSE 599b Project
# RNA Secondary Structure Design for Computational Aptamer Design

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-overview">Project Overview</a></li>
    <li><a href="#dependencies">Dependencies</a></li>
    <li><a href="#usage">Usage</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## Project Overview
The RNA secondary structure design problem aims to design an RNA sequence that folds into a target secondary structure with minimum free energy (MFE). This has major applications in nanotechnology with self-assembling structures, design of therapeutic agents (such as ribozymes), and will ultimately help to characterize RNA by function. To tackle this problem, this algorithm begins with an initialization step to generate random sequences for a particular structure and select the best initial sequence based on MFE and similarity to the target structure. Next, a stochastic local search algorithm is used to iteratively mutate the structure at the incorrectly folded bases. Finally, an improved sequence is returned.


## Dependencies
1. [NUPACK](http://www.nupack.org/)
2. [Pandas](https://anaconda.org/anaconda/pandas)
3. [Numpy](https://anaconda.org/anaconda/numpy)


<!-- USAGE EXAMPLES -->
## Usage

First, import the Python script:
  ```python
  from rna_inverse import *
  ```

### For a simple RNA secondary structure
First, input the target secondary structure in dot-parens notation and add primary constraints using the function `check_sequence_constraints`:
  ```python
target_structure = '(((((.....))))).....((......))......'
constraints = {}
constraints[1] = 'G'
constraints[8] = 'A'
# Convert constraints to full-list of constraints
constraints = check_sequence_constraints(target_structure, 
                                         constraints)
  ```

Next, generate the initial sequence using `initialize_sequences`:
  ```python
_, initial_sequence = initialize_sequences(target_structure, 50)
  ```

Evaluate the initial sequence using Nupack:
```python
nupack_structure, mfe = nupack_analyze_sequence(initial_sequence)
```

Mutate and improve the sequence using `mutate_sequence_iterate`. The function returns a dictionary detailing the success of designing a sequence, the number of iterations it took to reach an acceptable sequence, the final sequence, and the list of sequences the function moved through.
```python
number_iterations = 100
result = mutate_sequence_iterate(initial_sequence, 
                                 target_structure, 
                                 nupack_structure, 
                                 number_iterations)
improved_sequence = result['final_sequence']
```

