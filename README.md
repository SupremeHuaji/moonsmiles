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
test {
    // Validate SMILES string
    let is_valid = validate_smiles("C=C(C)O")  // true
    assert_eq(is_valid, true)
}
```

```moonbit
test {
    // Normalization
    let normalized_result = normalize_smiles("C = C ( C ) O")  // returns Result[String, String]
    let normalized = normalized_result.unwrap()  // or use match if .unwrap() is not available
    assert_eq(normalized, "C=C(C)O")
}
```

```moonbit
test {
    // Parse into molecular structure
    let molecule = parse_smiles("CCO")  // ethanol
    assert_eq(molecule.atoms.length(), 3)
}
```

```moonbit
test {
    // Generate canonical SMILES
    let molecule = parse_smiles("CCO")
    let canonical = generate_canonical_smiles(molecule)
    assert_eq(canonical, "CCO")
}
```

### Molecular Property Calculation
```moonbit
test {
    // Calculate molecular weight
    let molecule = parse_smiles("CCO")
    let mol_weight = calculate_molecular_weight(molecule)  // ~46.07 g/mol
    assert_eq(mol_weight > 46.0 && mol_weight < 46.1, true)
}
```

```moonbit
test {
    // Calculate LogP value (partition coefficient)
    let molecule = parse_smiles("CCO")
    let logp = calculate_logp(molecule)  // negative value indicates hydrophilicity
    assert_eq(logp < 0.0, true)
}
```

```moonbit
test {
    // Check drug-likeness
    let molecule = parse_smiles("CCO")
    let (passes_lipinski, _) = check_lipinski_rule_of_five(molecule)
    assert_eq(passes_lipinski, true)
}
```

```moonbit
test {
    // Count rotatable bonds
    let molecule = parse_smiles("CCO")
    let rotatable_bonds = count_rotatable_bonds(molecule)
    assert_eq(rotatable_bonds, 1)
}
```

### Structure Analysis
```moonbit
test {
    // Find all rings
    let molecule = parse_smiles("c1ccccc1")  // benzene
    let rings = find_all_rings(molecule)
    assert_eq(rings.length(), 1)
}
```

```moonbit
test {
    // Identify aromatic rings
    let molecule = parse_smiles("c1ccccc1")  // benzene
    let aromatic_rings = identify_aromatic_rings(molecule)
    assert_eq(aromatic_rings.length(), 1)
}
```

```moonbit
test {
    // Check for specific substructures
    let molecule = parse_smiles("CCO")
    let has_hydroxyl = contains_substructure(molecule, "O")  // check for hydroxyl group
    assert_eq(has_hydroxyl, true)
}
```

### Molecular Similarity
```moonbit
test {
    // Calculate similarity between two molecules
    let mol1 = parse_smiles("CCO")
    let mol2 = parse_smiles("CCC")
    let similarity = calculate_similarity(mol1, mol2)
    assert_eq(similarity > 0.4, true)
}
```

```moonbit
test {
    // Calculate similarity directly from SMILES
    let smiles_similarity = calculate_smiles_similarity("CCO", "CCC")
    assert_eq(smiles_similarity > 0.4, true)
}
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
