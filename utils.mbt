///|
/// Converts a string to lowercase
pub fn to_lower_case(s : String) -> String {
  let mut result = ""
  for c in s {
    let code = c.to_int()
    // ASCII uppercase letters are 65-90 (A-Z)
    if code >= 65 && code <= 90 {
      // 转为小写字母对应的字符
      let lowercase_char = Int::unsafe_to_char(code + 32).to_string()
      result = result + lowercase_char
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
      // 转为大写字母对应的字符
      let uppercase_char = Int::unsafe_to_char(code - 32).to_string()
      result = result + uppercase_char
    } else {
      result = result + c.to_string()
    }
  }
  return result
}

///|
/// Checks if a character is an uppercase letter (A-Z)
pub fn is_upper(c : Char) -> Bool {
  let code = c.to_int()
  // ASCII uppercase letters are 65-90 (A-Z)
  return code >= 65 && code <= 90
}

///|
/// 解析字符是否为小写字母
pub fn is_lower(c : Int) -> Bool {

  // ASCII lowercase letters are 97-122 (a-z)
  return c >= 97 && c <= 122
}

///|
/// 将数字字符串转换为整数
pub fn string_to_int(s : String) -> Int {
  let mut result = 0
  for c in s {
    let digit = c.to_int() - '0'.to_int()
    if digit >= 0 && digit <= 9 {
      result = result * 10 + digit
    } else {
      return 0 // 非数字字符，返回0
    }
  }
  return result
}

///|
pub fn is_lowercase(c : Int) -> Bool {
  c >= 97 && c <= 122 // ASCII 'a' 到 'z'
}
