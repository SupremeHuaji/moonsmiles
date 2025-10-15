//|

///|
/// Represents an atom in a molecule with enhanced properties
pub(all) struct Atom {
  symbol : String // 原子符号
  mut id : Int // 原子ID
  bonds : Array[Bond] // 连接的键
  charge : Int // 电荷
  chiral : String // 手性标记 (R/S)
  aromatic : Bool // 是否芳香
  isotope : Int? // 同位素质量数
  explicit_hydrogens : Int? // 显式氢原子数
} derive(Show)

///|
/// Represents a bond between two atoms
pub(all) struct Bond {
  atom1 : Int // ID of the first atom
  atom2 : Int // ID of the second atom
  mut bond_type : String
} derive(Show) // Bond type (e.g., "-", "=", "#")

///|
/// Represents a molecule
pub(all) struct Molecule {
  atoms : Array[Atom] // List of atoms in the molecule
  bonds : Array[Bond]
} derive(Show) // List of bonds in the molecule

///|
pub fn Molecule::new() -> Molecule {
  return Molecule::{
    atoms: Array::new(), // Initialize an empty array for atoms
    bonds: Array::new(),
  } // Initialize an empty array for bonds
}

///|
/// Creates a new Atom instance with default values
pub fn Atom::new(symbol : String, id : Int) -> Atom {
  {
    symbol,
    id,
    bonds: Array::new(),
    charge: 0,
    chiral: "",
    aromatic: is_lowercase(symbol[0]),
    isotope: None,
    explicit_hydrogens: None,
  }
}
