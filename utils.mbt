///|
/// Converts a string to lowercase
pub fn to_lower_case(s : String) -> String {
  let mut result = ""
  for c in s {
    let code = c.to_int()
    // ASCII uppercase letters are 65-90 (A-Z)
    if code >= 65 && code <= 90 {
      // Convert to lowercase by adding 32
      result = result + (code + 32).to_string()
    } else {
      result = result + c.to_string()
    }
  }
  return result
}

///|
/// Converts a string to uppercase
pub fn to_upper_case(s : String) -> String {
  let mut result = ""
  for c in s {
    let code = c.to_int()
    // ASCII lowercase letters are 97-122 (a-z)
    if code >= 97 && code <= 122 {
      // Convert to uppercase by subtracting 32
      result = result + (code - 32).to_string()
    } else {
      result = result + c.to_string()
    }
  }
  return result
}

///|
/// Checks if a character is a digit (0-9)
pub fn is_digit(c : Int) -> Bool {
  return c >= 48 && c <= 57
}

///|
/// Checks if a character is an uppercase letter (A-Z)
pub fn is_upper(c : String) -> Bool {
  if c.length() != 1 {
    return false
  }
  let code = c[0]
  // ASCII uppercase letters are 65-90 (A-Z)
  return code >= 65 && code <= 90
}

///|
/// Checks if a character is a lowercase letter (a-z)
pub fn is_lower(c : String) -> Bool {
  if c.length() != 1 {
    return false
  }
  let code = c[0]
  // ASCII lowercase letters are 97-122 (a-z)
  return code >= 97 && code <= 122
}
