///|
/// Test for validate_smiles
test "validate_smiles with valid SMILES" {
  assert_eq(validate_smiles("C=C(C)O"), true)
}

///|
/// Test for tokenize_smiles
test "tokenize_smiles splits SMILES into tokens" {
  let tokens = tokenize_smiles("C=C(C)O")
  assert_eq(tokens, ["C", "=", "C", "(", "C", ")", "O"])
}

///|
/// Test for count_atoms
test "count_atoms counts the number of atoms in SMILES" {
  assert_eq(count_atoms("C=C(C)O"), 4)
}

///|
/// Test for count_bonds
test "count_bonds counts the number of bonds in SMILES" {
  assert_eq(count_bonds("C=C(C)O"), 3)
}

///|
/// Test for has_balanced_parentheses
test "has_balanced_parentheses checks parentheses balance" {
  assert_eq(has_balanced_parentheses("C=C(C)O"), true)
  assert_eq(has_balanced_parentheses("C=C(C(O)"), false)
}

///|
/// Test for parse_smiles
test "parse_smiles parses SMILES into a Molecule" {
  let mol = parse_smiles("C=C(C)O")
  assert_eq(mol.atoms.length(), 4)
  assert_eq(mol.bonds.length(), 3)
}

///|
/// Test for generate_smiles
test "generate_smiles generates SMILES from a Molecule" {
  let mol = parse_smiles("C=C(C)O")
  let smiles = generate_smiles(mol)
  assert_eq(smiles, "C=C(C)O")
}

///|
/// Test for is_atom
test "is_atom identifies atomic symbols" {
  assert_eq(is_atom("C"), true)
  assert_eq(is_atom("O"), true)
  assert_eq(is_atom("N"), true)
  assert_eq(is_atom("H"), true)
  assert_eq(is_atom("X"), false) // 不是有效的原子符号
  assert_eq(is_atom("="), false) // 是键符号
}

///|
/// Test for is_bond
test "is_bond identifies bond symbols" {
  assert_eq(is_bond("-"), true)
  assert_eq(is_bond("="), true)
  assert_eq(is_bond("#"), true)
  assert_eq(is_bond("C"), false) // 是原子符号
  assert_eq(is_bond("("), false) // 是分支符号
}

///|
/// Test for is_branch
test "is_branch identifies branch symbols" {
  assert_eq(is_branch("("), true)
  assert_eq(is_branch(")"), true)
  assert_eq(is_branch("C"), false) // 是原子符号
  assert_eq(is_branch("-"), false) // 是键符号
}

///|
/// Test for remove_whitespace
test "remove_whitespace removes spaces from SMILES" {
  assert_eq(remove_whitespace("C = C ( C ) O"), "C=C(C)O")
  assert_eq(remove_whitespace("C\tC"), "CC")
}

///|
/// Test for normalize_smiles
test "normalize_smiles validates and cleans SMILES" {
  match normalize_smiles("C = C ( C ) O") {
    Ok(smiles) => assert_eq(smiles, "C=C(C)O")
    Err(_) => panic()
  }
  match normalize_smiles("C = C ( C") {
    Ok(_) => panic()
    Err(err) => assert_eq(err, "Unbalanced parentheses in SMILES string")
  }
}

///|
/// Test for validate_and_parse
test "validate_and_parse validates and parses SMILES" {
  match validate_and_parse("C=C(C)O") {
    Ok(mol) => {
      assert_eq(mol.atoms.length(), 4)
      assert_eq(mol.bonds.length(), 3)
    }
    Err(_) => panic()
  }
  match validate_and_parse("C=C(C") {
    Ok(_) => panic()
    Err(err) => assert_eq(err, "Unbalanced parentheses in SMILES string")
  }
}

///|
/// Test for count_atoms_in_molecule
test "count_atoms_in_molecule counts atoms in molecule" {
  let mol = parse_smiles("C=C(C)O")
  assert_eq(count_atoms_in_molecule(mol), 4)
}

///|
/// Test for count_bonds_in_molecule
test "count_bonds_in_molecule counts bonds in molecule" {
  let mol = parse_smiles("C=C(C)O")
  assert_eq(count_bonds_in_molecule(mol), 3)
}

///|
/// Test for is_ring_number
test "is_ring_number identifies ring numbers" {
  assert_eq(is_ring_number("1"), true)
  assert_eq(is_ring_number("9"), true)
  assert_eq(is_ring_number("0"), true)
  assert_eq(is_ring_number("C"), false)
  assert_eq(is_ring_number("-"), false)
}

///|
/// Test for tokenize_smiles_with_rings
test "tokenize_smiles_with_rings splits SMILES with ring numbers" {
  let tokens = tokenize_smiles_with_rings("C1=CC=CC=C1")
  assert_eq(tokens, ["C", "1", "=", "C", "C", "=", "C", "C", "=", "C", "1"])
}

///|
/// Test for calculate_molecular_weight
test "calculate_molecular_weight computes molecular weight" {
  let mol = parse_smiles("CCO") // 乙醇 C2H6O
  let weight = calculate_molecular_weight(mol)
  // C: 12.011 * 2 + H: 1.008 * 6 + O: 15.999 = 46.069
  assert_eq(weight > 46.0 && weight < 47.0, true)
}

///|
/// Test for calculate_smiles_weight
test "calculate_smiles_weight computes weight from SMILES" {
  let weight = calculate_smiles_weight("CCO") // 乙醇
  assert_eq(weight > 46.0 && weight < 47.0, true)
}

///|
/// Test for contains_any
test "contains_any checks if string contains any characters" {
  assert_eq(contains_any("C1=CC=CC=C1", ["1", "2"]), true)
  assert_eq(contains_any("C=CC=C", ["1", "2"]), false)
}

///|
/// Test for count_hydrogen_bond_acceptors
test "count_hydrogen_bond_acceptors counts H-bond acceptors" {
  // 甲醇 CH3OH 有1个氢键受体（氧原子）
  let methanol = parse_smiles("CO")
  assert_eq(count_hydrogen_bond_acceptors(methanol), 1)

  // 乙醚 CH3OCH3 有1个氢键受体
  let ether = parse_smiles("COC")
  assert_eq(count_hydrogen_bond_acceptors(ether), 1)
}

///|
/// Test for calculate_fingerprint
test "calculate_fingerprint generates molecular fingerprint" {
  let mol = parse_smiles("CCO")
  let fingerprint = calculate_fingerprint(mol)
  assert_eq(fingerprint.length(), 1024)

  // 至少有几位被设置为true
  let mut true_count = 0
  for bit in fingerprint {
    if bit {
      true_count += 1
    }
  }
  assert_eq(true_count > 0, true)
}

///|
/// Test for calculate_similarity
test "calculate_similarity computes Tanimoto similarity" {
  let mol1 = parse_smiles("CCO") // 乙醇
  let mol2 = parse_smiles("CCC") // 丙烷
  let sim = calculate_similarity(mol1, mol2)

  // 相似度应该在0到1之间
  assert_eq(sim >= 0.0 && sim <= 1.0, true)

  // 相同分子的相似度应为1.0
  let mol3 = parse_smiles("CCO")
  assert_eq(calculate_similarity(mol1, mol3), 1.0)
}

///|
/// Test for calculate_smiles_similarity
test "calculate_smiles_similarity computes similarity of SMILES" {
  let sim = calculate_smiles_similarity("CCO", "CCC")
  assert_eq(sim >= 0.0 && sim <= 1.0, true)
  let self_sim = calculate_smiles_similarity("CCO", "CCO")
  assert_eq(self_sim, 1.0)
}
