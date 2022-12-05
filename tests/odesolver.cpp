#include <iostream>
#include <gtest/gtest.h>
#include "../include/odesolver.h"
#include "../include/utility.h"

std::pair<parser::Lexeme, std::string> pair(parser::Lexeme lexeme, std::string string) {
    return std::make_pair(lexeme, std::move(string));
}

TEST(ParserTest, Basic) {
    using namespace parser;
    std::stringstream ss("x' = x");

    Trie trie;
    trie.addString("x");

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(ParserTest, Numbers) {
    using namespace parser;
    std::stringstream ss("x' = 5 x");

    Trie trie;
    trie.addString("x");

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "5"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(ParserTest, Complex) {
    using namespace parser;
    std::stringstream ss("x' = - 7 a y + 3 z");

    Trie trie;
    for (const auto& string: {"x", "y", "z", "a"}) {
        trie.addString(string);
    }

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Minus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "7"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "a"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "y"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Plus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "3"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "z"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
