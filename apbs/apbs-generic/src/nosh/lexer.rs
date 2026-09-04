// APBS nosh/lexer - Tokenizer for input file

/// Token types
#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    Keyword(String),
    String(String),
    Integer(i32),
    Real(f64),
    Lparen,
    Rparen,
    Semicolon,
    Equal,
    Comma,
    Eof,
}

/// Tokenizer
pub struct Lexer {
    input: Vec<char>,
    pos: usize,
}

impl Lexer {
    pub fn new(input: &str) -> Self {
        Self {
            input: input.chars().collect(),
            pos: 0,
        }
    }

    /// Get next token
    pub fn next_token(&mut self) -> Token {
        self.skip_whitespace();

        if self.pos >= self.input.len() {
            return Token::Eof;
        }

        match self.input[self.pos] {
            '(' => { self.pos += 1; Token::Lparen }
            ')' => { self.pos += 1; Token::Rparen }
            ';' => { self.pos += 1; Token::Semicolon }
            '=' => { self.pos += 1; Token::Equal }
            ',' => { self.pos += 1; Token::Comma }
            '"' | '\'' => self.read_string(),
            c if c.is_ascii_digit() || c == '-' || c == '+' => self.read_number(),
            c if c.is_alphabetic() || c == '_' => self.read_word(),
            _ => {
                self.pos += 1;
                Token::Eof
            }
        }
    }

    fn skip_whitespace(&mut self) {
        while self.pos < self.input.len() && self.input[self.pos].is_whitespace() {
            self.pos += 1;
        }
    }

    fn read_string(&mut self) -> Token {
        let quote = self.input[self.pos];
        self.pos += 1;
        let start = self.pos;
        while self.pos < self.input.len() && self.input[self.pos] != quote {
            self.pos += 1;
        }
        let s: String = self.input[start..self.pos].iter().collect();
        if self.pos < self.input.len() {
            self.pos += 1; // skip closing quote
        }
        Token::String(s)
    }

    fn read_number(&mut self) -> Token {
        let start = self.pos;
        if self.input[self.pos] == '-' || self.input[self.pos] == '+' {
            self.pos += 1;
        }
        while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos < self.input.len() && self.input[self.pos] == '.' {
            self.pos += 1;
            while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
            let s: String = self.input[start..self.pos].iter().collect();
            Token::Real(s.parse().unwrap_or(0.0))
        } else {
            let s: String = self.input[start..self.pos].iter().collect();
            Token::Integer(s.parse().unwrap_or(0))
        }
    }

    fn read_word(&mut self) -> Token {
        let start = self.pos;
        while self.pos < self.input.len()
            && (self.input[self.pos].is_alphanumeric() || self.input[self.pos] == '_')
        {
            self.pos += 1;
        }
        let s: String = self.input[start..self.pos].iter().collect();
        Token::Keyword(s)
    }
}
