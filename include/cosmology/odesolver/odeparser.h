#pragma once

#include <istream>
#include <unordered_set>
#include <unordered_map>
#include <string>

namespace parser {
    class Trie {
    private:
        struct Node {
            std::unordered_map<char, Node *> children;

            ~Node() {
                for (const auto &[ch, child]: children) {
                    delete child;
                }
            }
        };

    public:
        ~Trie() {
            delete root_;
        }

        void addString(const std::string &string) {
            Node *node = root_;

            for (char c: string) {
                if (!node->children.count(c)) {
                    node->children[c] = new Node;
                }
                node = node->children[c];
            }
        }

        std::string consumeStream(std::istream &stream) const {
            Node *node = root_;
            std::string string;

            while (stream) {
                char c = static_cast<char>(stream.peek());

                if (node->children.count(c)) {
                    string.push_back(c);
                    node = node->children[c];
                    stream >> c;
                } else {
                    break;
                }
            }

            return string;
        }

    private:
        Node *root_ = new Node;
    };

    enum Lexeme {
        Left, Right, Plus, Minus, Derivative, LeftBracket, RightBracket, Number, End
    };

    class EquationLexer {
    public:
        EquationLexer(const Trie &trie, std::istream &stream) : trie_(trie), stream_(stream) {
        }

        std::pair<Lexeme, std::string> next() {
            while (stream_) {
                char c = static_cast<char>(stream_.peek());
                if (c == EOF) {
                    break;
                }
                if (specialSymbols_.count(c)) {
                    stream_.get(c);
                    if (c == '+') {
                        return {Lexeme::Plus, ""};
                    } else if (c == '-') {
                        return {Lexeme::Minus, ""};
                    } else if (c == '(') {
                        return {Lexeme::LeftBracket, ""};
                    } else if (c == ')') {
                        return {Lexeme::RightBracket, ""};
                    } else if (c == '\'') {
                        if (!left_) {
                            throw std::logic_error("can't have derivative in the right part of an equation ");
                        }
                        return {Lexeme::Derivative, ""};
                    } else if (c == '=') {
                        if (!left_) {
                            throw std::logic_error("can't have more than one '=' in an equation");
                        }
                        left_ = false;
                    }
                } else if (c >= '0' && c <= '9') {
                    return {Lexeme::Number, consumeNumber()};
                } else {
                    auto word = trie_.consumeStream(stream_);
                    if (left_) {
                        return {Lexeme::Left, word};
                    } else {
                        return {Lexeme::Right, word};
                    }
                }
            }
            return {Lexeme::End, ""};
        }

    private:
        std::string consumeNumber() {
            std::string number;

            while (stream_) {
                char c = static_cast<char>(stream_.peek());

                if (c >= '0' && c <= '9') {
                    stream_ >> c;
                    number.push_back(c);
                } else {
                    break;
                }
            }

            return number;
        }

    private:
        std::istream &stream_;
        const Trie &trie_;
        inline static const std::unordered_set<char> specialSymbols_{'+', '-', '*', '=', '\'', ' ', '(', ')'};
        bool left_ = true;
    };
} // namespace parser
