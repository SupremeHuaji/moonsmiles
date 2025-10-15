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
  assert_eq(fingerprint.length(), 2048)

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

///|
/// Test for generate_formula
test "generate_formula creates molecular formula" {
  // 乙醇测试
  let ethanol = parse_smiles("CCO")
  assert_eq(generate_formula(ethanol), "C2H6O")

  // 苯测试
  let benzene = parse_smiles("c1ccccc1")
  assert_eq(generate_formula(benzene), "C6H6")

  // 苯酚测试
  let phenol = parse_smiles("c1ccccc1O")
  assert_eq(generate_formula(phenol), "C6H6O")
}

///|
/// Test for count_hydrogen_bond_donors
test "count_hydrogen_bond_donors counts H-bond donors" {
  // 甲醇 CH3OH 有1个氢键供体
  let methanol = parse_smiles("CO")
  assert_eq(count_hydrogen_bond_donors(methanol), 1)

  // 乙二醇有两个氢键供体
  let glycol = parse_smiles("OCCO")
  assert_eq(count_hydrogen_bond_donors(glycol), 2)

  // 苯 C6H6 没有氢键供体
  let benzene = parse_smiles("c1ccccc1")
  assert_eq(count_hydrogen_bond_donors(benzene), 0)
}

///|
/// Test for contains_substructure
test "contains_substructure checks for substructure patterns" {
  // 乙醇含有甲基
  let ethanol = parse_smiles("CCO")
  assert_eq(contains_substructure(ethanol, "C"), true)

  // 乙醇含有羟基
  assert_eq(contains_substructure(ethanol, "O"), true)

  // 乙醇不含苯环
  assert_eq(contains_substructure(ethanol, "c1ccccc1"), false)

  // 检查酯基
  let methyl_acetate = parse_smiles("CC(=O)OC")
  assert_eq(contains_substructure(methyl_acetate, "C(=O)O"), true)
}

///|
/// Test for calculate_logp
test "calculate_logp computes LogP value" {
  // 测试乙醇
  let ethanol = parse_smiles("CCO")
  let ethanol_logp = calculate_logp(ethanol)

  assert_eq(ethanol_logp < 0.0, true) // 乙醇LogP应为负值(亲水性)

  // 测试正辛烷
  let octane = parse_smiles("CCCCCCCC")
  let octane_logp = calculate_logp(octane)
  assert_eq(octane_logp > 3.0, true) // 正辛烷LogP应为正值(疏水性)
}

///|
/// Test for check_lipinski_rule_of_five
test "check_lipinski_rule_of_five evaluates druglikeness" {
  // 测试阿司匹林 (符合Lipinski规则)
  let aspirin = parse_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
  let (pass, _) = check_lipinski_rule_of_five(aspirin)
  assert_eq(pass, true)

  // 测试环孢素 (不符合Lipinski规则 - 分子量大)
  // 简化版环孢素模型，实际结构更复杂
  let cyclosporine = parse_smiles(
    "CC1C(C(C(N(C(=O)C(N(C(=O)C(N(C(=O)C(N(C(=O)C(N(C1=O)C(C(C)C)C)C)C)C)C)C)CC2=CC=CC=C2)C)C)C)C",
  )
  let (pass2, violations) = check_lipinski_rule_of_five(cyclosporine)
  assert_eq(pass2, false)
  assert_eq(violations.length() > 0, true)
}

///|
/// Test for tokenize_smiles_extended
test "tokenize_smiles_extended handles bracketed atoms" {
  let tokens = tokenize_smiles_extended("C[CH3][OH+]")
  assert_eq(tokens.length(), 3)
  assert_eq(tokens[0], "C")
  assert_eq(tokens[1], "[CH3]")
  assert_eq(tokens[2], "[OH+]")
}

///|
/// Test for parse_bracketed_atom
test "parse_bracketed_atom parses atom properties" {
  let oh_plus = parse_bracketed_atom("[OH+]")
  assert_eq(oh_plus.symbol, "O")
  assert_eq(oh_plus.charge, 1)
  let c13 = parse_bracketed_atom("[13CH3]")
  assert_eq(c13.symbol, "C")
  assert_eq(c13.isotope, Some(13))
}

///|
/// Test for get_valence and get_normal_valence
test "get_valence returns correct valence values" {
  // 创建测试原子
  let carbon = Atom::new("C", 0)
  let nitrogen = Atom::new("N", 1)
  let oxygen = Atom::new("O", 2)
  let aromatic_c = Atom::new("c", 3)
  assert_eq(get_normal_valence("C"), 4)
  assert_eq(get_normal_valence("N"), 3)
  assert_eq(get_normal_valence("O"), 2)
  assert_eq(get_valence(carbon), 4)
  assert_eq(get_valence(nitrogen), 3)
  assert_eq(get_valence(oxygen), 2)
  assert_eq(get_valence(aromatic_c), 3) // 芳香碳价态为3
}

///|
/// Test for count_rotatable_bonds
test "count_rotatable_bonds counts rotatable bonds" {
  // 丙醇 CH3CH2CH2OH 有2个可旋转键
  let propanol = parse_smiles("CCCO")
  assert_eq(count_rotatable_bonds(propanol), 2)

  // 二乙基醚 CH3CH2OCH2CH3 有2个可旋转键
  let diethyl_ether = parse_smiles("CCOCC")
  assert_eq(count_rotatable_bonds(diethyl_ether), 2)

  // 苯 C6H6 没有可旋转键
  let benzene = parse_smiles("c1ccccc1")
  assert_eq(count_rotatable_bonds(benzene), 0)
}

///|
/// Test for find_all_rings
test "find_all_rings identifies all rings in molecule" {
  // 苯只有一个环
  let benzene = parse_smiles("c1ccccc1")
  let benzene_rings = find_all_rings(benzene)
  assert_eq(benzene_rings.length(), 1)
  assert_eq(benzene_rings[0].length(), 6)

  // 萘有两个环
  let naphthalene = parse_smiles("c1ccc2ccccc2c1")
  let naphthalene_rings = find_all_rings(naphthalene)
  assert_eq(naphthalene_rings.length() >= 2, true)
}

///|
/// Test for identify_aromatic_rings
test "identify_aromatic_rings finds aromatic rings" {
  // 苯是芳香环
  let benzene = parse_smiles("c1ccccc1")
  let aromatic_rings = identify_aromatic_rings(benzene)
  assert_eq(aromatic_rings.length(), 1)

  // 环己烷不是芳香环
  let cyclohexane = parse_smiles("C1CCCCC1")
  let non_aromatic = identify_aromatic_rings(cyclohexane)
  assert_eq(non_aromatic.length(), 0)
}

///|
/// Test for calculate_implicit_hydrogens
test "calculate_implicit_hydrogens counts implicit hydrogens" {
  // 甲烷 (CH4)
  let methane_mol = parse_smiles("C")
  let methane_atom = methane_mol.atoms[0]
  assert_eq(calculate_implicit_hydrogens(methane_atom, methane_mol), 4)

  // 乙烯 (C2H4) - 碳参与双键
  let ethylene_mol = parse_smiles("C=C")
  let ethylene_c1 = ethylene_mol.atoms[0]
  assert_eq(calculate_implicit_hydrogens(ethylene_c1, ethylene_mol), 2)
}

///|
/// Test for generate_canonical_smiles
test "generate_canonical_smiles creates standardized SMILES" {
  // 不同表示的乙醇应该生成相同的规范形式
  let ethanol1 = parse_smiles("CCO")
  let ethanol2 = parse_smiles("OCC")
  let canonical1 = generate_canonical_smiles(ethanol1)
  let canonical2 = generate_canonical_smiles(ethanol2)
  assert_eq(canonical1, canonical2)

  // 丙醇不同表示
  let propanol1 = parse_smiles("CCCO")
  let propanol2 = parse_smiles("OCCC")
  assert_eq(
    generate_canonical_smiles(propanol1),
    generate_canonical_smiles(propanol2),
  )
}

///|
/// Test for is_lowercase and to_lower_case functions
test "string utility functions work correctly" {
  assert_eq(to_lower_case("SMILES"), "smiles")
  assert_eq(to_upper_case("smiles"), "SMILES")
}

///|
/// Test for string_to_int function
test "string_to_int converts strings to integers" {
  assert_eq(string_to_int("123"), 123)
  assert_eq(string_to_int("0"), 0)
  assert_eq(string_to_int("42"), 42)
}

///|
/// Test for ring detection in complex molecules
test "ring detection in complex structures" {
  // 测试联苯
  let biphenyl = parse_smiles_with_rings("c1ccccc1-c2ccccc2")
  let biphenyl_rings = find_all_rings(biphenyl)
  assert_eq(biphenyl_rings.length(), 2)

  // 测试环己烷
  let cyclohexane = parse_smiles_with_rings("C1CCCCC1")
  let ring_bonds = identify_ring_bonds(cyclohexane)
  assert_eq(ring_bonds.size(), 6) // 环己烷有6个环键
}

///|
/// Test for atom property handling
test "atom properties are correctly processed" {
  // 测试带电荷的原子
  let ammonium = parse_smiles("[NH4+]")
  let nitrogen = ammonium.atoms[0]
  assert_eq(nitrogen.charge, 1)

  // 测试隐式氢
  let water = parse_smiles("O")
  let oxygen = water.atoms[0]
  assert_eq(calculate_implicit_hydrogens(oxygen, water), 2)
}
