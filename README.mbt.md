# MoonSMILES - A Molecular Processing Library for MoonBit

MoonSMILES is a cheminformatics library implemented in MoonBit language, focusing on processing SMILES (Simplified Molecular Input Line Entry System) strings and calculating molecular properties. The library provides a suite of tools for molecular structure analysis, property prediction, and similarity comparison.

## Key Features

### SMILES String Processing
- **Validation & Parsing**: Validate SMILES string correctness, check bracket balance, and parse into internal molecular representation
- **Normalization & Standardization**: Clean SMILES strings and generate canonical forms to ensure consistent representation
- **Format Conversion**: Convert between SMILES strings and molecular structure objects

### Molecular Structure Analysis
- **Atom & Bond Counting**: Calculate the number of atoms and bonds in molecules
- **Ring Structure Detection**: Find all ring structures and distinguish between aromatic and non-aromatic rings
- **Substructure Searching**: Detect specific structural patterns within molecules
- **Formula Generation**: Generate molecular formulas (like C2H6O) from structure

### Drug Property Analysis
- **Molecular Weight Calculation**: Compute exact mass of molecules
- **LogP Prediction**: Estimate octanol-water partition coefficient (lipophilicity/hydrophilicity indicator)
- **Lipinski Rule Checking**: Evaluate drug-likeness of molecules
- **Hydrogen Bond Properties**: Calculate hydrogen bond donors and acceptors
- **Rotatable Bond Counting**: Analyze molecular flexibility

### Molecular Similarity Comparison
- **Molecular Fingerprint Generation**: Create bit vectors representing molecular structural features
- **Similarity Calculation**: Compare structural similarity between molecules using Tanimoto coefficient

## Usage Examples

### Basic SMILES Processing
```moonbit
    // Validate SMILES string
    let is_valid = validate_smiles("C=C(C)O")  // true

    // Normalization
    let normalized = normalize_smiles("C = C ( C ) O")  // returns "C=C(C)O"

    // Parse into molecular structure
    let molecule = parse_smiles("CCO")  // ethanol

    // Generate canonical SMILES
    let canonical = generate_canonical_smiles(molecule)
```

### Molecular Property Calculation
```moonbit
    // Calculate molecular weight
    let mol_weight = calculate_molecular_weight(molecule)  // ~46.07 g/mol

    // Calculate LogP value (partition coefficient)
    let logp = calculate_logp(molecule)  // negative value indicates hydrophilicity

    // Check drug-likeness
    let (passes_lipinski, violations) = check_lipinski_rule_of_five(molecule)

    // Count rotatable bonds
    let rotatable_bonds = count_rotatable_bonds(molecule)
```

### Structure Analysis
```moonbit
    // Find all rings
    let rings = find_all_rings(parse_smiles("c1ccccc1"))  // benzene ring

    // Identify aromatic rings
    let aromatic_rings = identify_aromatic_rings(molecule)

    // Check for specific substructures
    let has_hydroxyl = contains_substructure(molecule, "O")  // check for hydroxyl group
```

### Molecular Similarity
```moonbit
    // Calculate similarity between two molecules
    let similarity = calculate_similarity(mol1, mol2)

    // Calculate similarity directly from SMILES
    let smiles_similarity = calculate_smiles_similarity("CCO", "CCC")
```

## Special Features

- Support for extended atom representations including isotope labeling and charge states: `[13CH3]`, `[NH4+]`
- Support for aromaticity representation: `c1ccccc1` (benzene)
- Handling of complex ring structures and multi-ring systems

## Application Domains

- Drug discovery and design
- Compound screening and virtual screening
- Chemical database searching
- Structure-activity relationship (SAR) studies
- Chemical education and research

MoonSMILES library provides a powerful and easy-to-use toolset for cheminformatics and computational chemistry research, implemented in MoonBit language for efficiency, safety, and readability.