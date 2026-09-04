// APBS vstring.rs - String utilities
// Port of src/generic/vstring.h / vstring.c

/// Case-insensitive string comparison (like strcasecmp)
pub fn strcasecmp(s1: &str, s2: &str) -> std::cmp::Ordering {
    s1.to_ascii_lowercase().cmp(&s2.to_ascii_lowercase())
}

/// Case-insensitive string equality
pub fn strcasecmp_eq(s1: &str, s2: &str) -> bool {
    s1.eq_ignore_ascii_case(s2)
}

/// Test if entire string consists of ASCII digits
pub fn isdigit(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// Word-wrapping text formatter.
/// Wraps text to `right_margin` columns with `left_padding` spaces.
pub fn wrappedtext(str: &str, right_margin: usize, left_padding: usize) -> String {
    let mut result = String::new();
    let padding: String = " ".repeat(left_padding);
    let mut col = 0;

    for word in str.split_whitespace() {
        if col + word.len() + 1 > right_margin && col > 0 {
            result.push('\n');
            result.push_str(&padding);
            col = 0;
        }
        if col > 0 {
            result.push(' ');
            col += 1;
        }
        result.push_str(word);
        col += word.len();
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strcasecmp() {
        assert_eq!(strcasecmp("abc", "ABC"), std::cmp::Ordering::Equal);
        assert_eq!(strcasecmp("abc", "abd"), std::cmp::Ordering::Less);
        assert_eq!(strcasecmp("abd", "abc"), std::cmp::Ordering::Greater);
    }

    #[test]
    fn test_isdigit() {
        assert!(isdigit("12345"));
        assert!(!isdigit("123a5"));
        assert!(!isdigit(""));
        assert!(!isdigit("abc"));
    }
}
