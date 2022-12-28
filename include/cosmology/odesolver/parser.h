#pragma once

#include <istream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <string>
#include <utility>

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
        Variable, Number,
        Plus, Minus, Multiplication, Division, Derivative,
        Equal, LeftBracket, RightBracket, End
    };

    using Token = std::pair<Lexeme, std::string>;

    class RawExpressionLexer {
    public:
        RawExpressionLexer(const Trie &trie, std::istream &stream) : trie_(trie), stream_(stream) {
        }

        Token next() {
            while (stream_) {
                char c = static_cast<char>(stream_.peek());
                if (c == EOF) {
                    break;
                }
                if (specialSymbols_.count(c)) {
                    stream_.get(c);
                    switch (c) {
                        case '+':
                            return {Lexeme::Plus, ""};
                        case '-':
                            return {Lexeme::Minus, ""};
                        case '*':
                            return {Lexeme::Multiplication, ""};
                        case '/':
                            return {Lexeme::Division, ""};
                        case '(':
                            return {Lexeme::LeftBracket, ""};
                        case ')':
                            return {Lexeme::RightBracket, ""};
                        case '\'':
                            return {Lexeme::Derivative, ""};
                        case '=':
                            return {Lexeme::Equal, ""};
                        default:
                            continue;
                    }
                } else if (c >= '0' && c <= '9') {
                    return {Lexeme::Number, consumeNumber()};
                } else {
                    return {Lexeme::Variable, trie_.consumeStream(stream_)};
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
        inline static const std::unordered_set<char> specialSymbols_{'+', '-', '*', '/', '=', '\'', ' ', '(', ')'};
    };

    class ExpressionLexer {
    public:
        ExpressionLexer(const Trie &trie, std::istream &stream) : lexer_(trie, stream) {
        }

        Token next() {
            if (valueQueue_.empty()) {
                Token token = lexer_.next();
                if (!isBatchStart(token)) {
                    return token;
                } else {
                    valueQueue_.push(token);
                }
            }
            if (valueQueue_.size() == 1 && isBatchStart(valueQueue_.back())) {
                while (true) {
                    Token token = lexer_.next();
                    if (isBatchContinuation(token)) {
                        valueQueue_.push({Lexeme::Multiplication, ""});
                        valueQueue_.push(token);
                    } else {
                        valueQueue_.push(token);
                        break;
                    }
                }
            }
            if (!valueQueue_.empty()) {
                Token token = valueQueue_.front();
                valueQueue_.pop();
                return token;
            } else {
                return {Lexeme::End, ""};
            }
        }

    private:
        static bool isValue(const Token &token) {
            return token.first == Lexeme::Variable || token.first == Lexeme::Number;
        }

        static bool isBatchStart(const Token &token) {
            return isValue(token) || token.first == Lexeme::RightBracket;
        }

        bool isBatchContinuation(const Token &token) {
            if (isValue(valueQueue_.back()) && (isValue(token) || token.first == Lexeme::LeftBracket)) {
                return true;
            }
            if (valueQueue_.back().first == Lexeme::RightBracket && isValue(token)) {
                return true;
            }
            return false;
        }

    private:
        RawExpressionLexer lexer_;
        std::queue<Token> valueQueue_;
    };

    class SyntaxTree {
    private:
        struct Node {
            Token token;
            int level;
            Node *left;
            Node *right;

            template<class... Args>
            explicit Node(int new_level, Args &&... args) : token(Token(std::forward<Args>(args)...)), level(new_level),
                                                            left(nullptr), right(nullptr) {
            }

            ~Node() {
                delete left;
                delete right;
            }
        };

    public:
        using Node = Node;

        ~SyntaxTree() {
            delete root_;
        }

        Node *getRoot() const {
            return root_;
        }

        Node *releaseRoot() {
            Node *root = root_;
            root_ = nullptr;
            return root;
        }

        void build(ExpressionLexer &lexer) {
            buildUntilAt(lexer, Lexeme::End, 0);
        }

        friend std::vector<Node*> getVertexNeighbors(const SyntaxTree& /*tree*/, Node* node) {
            std::vector<Node*> neighbors;

            for (auto neighbor: {node->right, node->left}) {
                if (neighbor) {
                    neighbors.push_back(neighbor);
                }
            }

            return neighbors;
        }

    private:
        void buildUntilAt(ExpressionLexer &lexer, Lexeme lexeme, int level) {
            Token token = lexer.next();

            while (token.first != lexeme) {
                if (token.first == Lexeme::LeftBracket) {
                    SyntaxTree bracketTree;
                    bracketTree.buildUntilAt(lexer, Lexeme::RightBracket, level + 1);
                    addNode(bracketTree.releaseRoot());
                } else {
                    Node *node = new Node(level, token);
                    addNode(node);
                }
                token = lexer.next();
            }
        }

        void addNode(Node *node) {
            Node *current = root_;
            Node *parent = nullptr;

            while (less(node, current)) {
                parent = current;
                current = parent->right;
            }

            if (current) {
                node->left = current;
            }
            if (!parent) {
                root_ = node;
            } else {
                parent->right = node;
            }
        }

        static int order(Node *node) {
            switch (node->token.first) {
                case Lexeme::Variable:
                case Lexeme::Number:
                    return 0;
                case Lexeme::Multiplication:
                case Lexeme::Division:
                    return 2;
                case Lexeme::Plus:
                    return 3;
                case Lexeme::Minus: {
                    if (node->left) {
                        return 3;
                    }
                    return 1;
                }
                default:
                    throw std::logic_error("unsupported lexeme in expression");
            }
        }

        static bool less(Node *lhs, Node *rhs) {
            if (!rhs) {
                return false;
            }
            if (!lhs) {
                return true;
            }
            return lhs->level > rhs->level || (lhs->level == rhs->level && order(lhs) < order(rhs));
        }

    private:
        Node *root_ = nullptr;
    };
} // namespace parser
