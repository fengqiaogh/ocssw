
#ifndef FOCS_PRODUCT
#define FOCS_PRODUCT


#include "StringUtils.hpp"

#include <iostream>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>
#include <memory>

// #include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/variant_fwd.hpp>

namespace focs {
    // TODO: change attributes to use std::multimap, which will take care of the vector crap
    // (assuming I don't need more than one dimension for an attribute)?
    // For now, the NetCdf reader just makes a comma-separated string out of non-scalar attributes.
    // Stuff downstream can use the StringUtils::stov functions to undo it.


    using AttributeType = boost::variant<char, std::string, float, double, long double,
        int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t>;
        // std::vector<std::string>, std::vector<float>, std::vector<double>, std::vector<long double>,
        // std::vector<int8_t>, std::vector<int16_t>, std::vector<int32_t>, std::vector<int64_t>,
        // std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>, std::vector<uint64_t>>;


    class Attribute {
        public:
            Attribute(const std::string& name, AttributeType value) : name_{name}, value_{value}{}

            const std::string& name() const {return name_;}
            const AttributeType& value() const {return value_;}
            const AttributeType& operator*() const {return value_;}

            double as_double() const {
                return boost::apply_visitor(AsDouble{}, value_);
            }

            std::string as_string() const {
                return boost::apply_visitor(AsString{}, value_);
            }

            template<typename T>
            std::vector<T> as_vector() const {
                return boost::apply_visitor(AsVector<T>{}, value_);
            }

            friend std::ostream& operator<<(std::ostream& out, const Attribute& attribute){
                out << attribute.name_ << "=";
                boost::apply_visitor(Print{out}, attribute.value_);
                return out;
            }

            friend bool operator<(const Attribute& left, const Attribute& right){
                // if (left.name_ != right.name_){
                if (!boost::iequals(left.name_, right.name_)){
                    return left.name_ < right.name_;
                }
                return boost::apply_visitor(LessThan{}, left.value_, right.value_);
            }

            friend bool operator==(const Attribute& left, const std::string& name){
                // return (left.name_ == name);
                return boost::iequals(left.name_, name);
            }
            friend bool operator==(const Attribute& left, const Attribute& right){
                // if (left.name_ != right.name_){
                if (!boost::iequals(left.name_, right.name_)){
                    return false;
                }
                return boost::apply_visitor(Equals{}, left.value_, right.value_);
            }
            static Attribute parse(const std::string& input){
                auto equals = input.find_first_of('=');
                size_t name_start = 0;
                while (input[name_start] == ' '){
                    name_start++;
                }
                size_t name_end = equals;
                while (input[name_end - 1] == ' '){
                    name_end--;
                }
                std::string name(input, name_start, name_end - name_start);

                while (input[++equals] == ' ');

                switch (input[equals]){
                    case '"':
                        return Attribute{name, input.substr(equals+1, input.size() - equals - 2)};
                    case '\'':
                        return Attribute{name, input[equals+1]};
                    default:
                        // TODO
                        // This should probably utilize value suffixes, like C (e.g., to force float, 5f)
                        // I also should check for bad parses
                        if (input.substr(equals) == "true"){
                            return Attribute{name, true};
                        } else if (input.substr(equals) == "false"){
                            return Attribute{name, false};
                        } else if (input.find_first_of('.') == std::string::npos){
                            return Attribute{name, std::stoi(input.substr(equals))};
                        } else {
                            return Attribute{name, std::stof(input.substr(equals))};
                        }
                }
            }


        private:
            std::string name_;
            AttributeType value_;
            struct Print : public boost::static_visitor<void> {
                public:
                    Print(std::ostream& out_) : out{out_} {}

                    void operator()(std::string s) const {
                        out << '"' << s << '"';
                    }
                    void operator()(uint8_t i) const {
                        out << static_cast<int>(i);
                    }
                    void operator()(int8_t i) const {
                        out << static_cast<int>(i);
                    }
                    void operator()(char c) const {
                        out << '\'' << c << '\'';
                    }
                    template <typename T>
                    void operator()(T s) const {
                        out << s;
                    }
                private:
                    std::ostream& out{std::cout};
            };
            struct AsDouble : public boost::static_visitor<double> {
                public:
                    AsDouble() {}
                    template <typename T>
                    double operator()(T s) const {
                        return s;
                    }
                    double operator()(std::string s) const {
                        return std::stod(s);
                    }
            };
            struct AsString : public boost::static_visitor<std::string> {
                public:
                    AsString() {}
                    template <typename T>
                    std::string operator()(T s) const {
                        return std::to_string(s);
                    }
                    std::string operator()(std::string s) const {
                        return s;
                    }
            };
            template<typename T>
            struct AsVector : public boost::static_visitor<std::vector<T>> {
                public:
                    AsVector() {}
                    template <typename I>
                    std::vector<T> operator()(I s) const {
                        return std::vector<T>(s);
                    }
                    std::vector<T> operator()(std::string s) const {
                        return focs::StringUtils::stov<T>(s);
                    }
            };

            struct Equals : public boost::static_visitor<bool> {
                template <typename T>
                bool operator() (T a, T b) const { return a == b; }
                bool operator()(const std::string& a, const std::string& b) const { return a == b; }

                template <typename T>
                bool operator()(std::string, T) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }
                template <typename T>
                bool operator()(T, std::string) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }

                // Integer comparison (extra work to avoid signed/unsigned comparison)
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value == std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a == b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && !std::is_signed<T>::value && std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (b < 0){
                        return false;
                    }
                    return static_cast<U>(a) == b;
                }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value && !std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (a < 0){
                        return false;
                    }
                    return a == static_cast<T>(b);
                }

                // Non-integer comparison
                template <typename T, typename U, typename std::enable_if<std::is_floating_point<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(a - b) < 0.000001; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(static_cast<double>(a) - b) < 0.000001; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<U>::value && std::is_floating_point<T>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(a - static_cast<double>(b)) < 0.000001; }
            };
            struct LessThan : public boost::static_visitor<bool> {
                template <typename T>
                bool operator() (T a, T b) const { return a < b; }
                bool operator()(const std::string& a, const std::string& b) const { return a < b; }

                template <typename T>
                bool operator()(std::string, T) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }
                template <typename T>
                bool operator()(T, std::string) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }

                // Integer comparison (extra work to avoid signed/unsigned comparison)
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value == std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && !std::is_signed<T>::value && std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (b < 0){
                        return false;
                    }
                    return static_cast<U>(a) < b;
                }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value && !std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (a < 0){
                        return true;
                    }
                    return a < static_cast<T>(b);
                }

                // Non-integer comparison
                template <typename T, typename U, typename std::enable_if<std::is_floating_point<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return static_cast<double>(a) < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<U>::value && std::is_floating_point<T>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < static_cast<double>(b); }
            };
    };
    class AttributeCondition {
        public:
            virtual ~AttributeCondition(){}
            virtual bool matches(const Attribute& attribute) const {(void)attribute; return false;}
            friend std::ostream& operator<<(std::ostream& out, const AttributeCondition&){return out << "unknown condition";}
    };
    class AttributeWild : public AttributeCondition {
        public:
            AttributeWild(const std::string& name) : name_{name} {}
            virtual ~AttributeWild() override {}
            bool matches(const Attribute& attribute) const override {
                return attribute.name() == name_;
            }
            friend std::ostream& operator<<(std::ostream& out, const AttributeWild& condition){
                return out << condition.name_ << "=*";
            }
        private:
            std::string name_;
    };
    template<typename T=double>
    class AttributeRange : public AttributeCondition {
        public:
            // AttributeRange(const std::string& name, T min, T max) : name_{name}, min_{min}, max_{max} {}
            AttributeRange(const std::string& name, T min, T max) : name_{name}, between_{Between<T>(min, max)} {}
            virtual ~AttributeRange() override {}
            bool matches(const Attribute& attribute) const override {
                return attribute.name() == name_ && boost::apply_visitor(between_, attribute.value());
            }
            // friend std::ostream& operator<<(std::ostream& out, const AttributeRange& condition){
            //     // return out name_ << ": " << condition.between_.min_ << " - " << condition.between_.max_;
            //     return out name_ << " between stuff";
            // }
        private:
            template<typename B>
            struct Between : public boost::static_visitor<bool> {
                Between(B min, B max) : min_{min}, max_{max} {}
                B min_;
                B max_;

                bool operator()(const std::string&) const { return false; }
                template <typename V>
                bool operator() (V a) const { return !(a < min_ || a > max_); }
            };

            std::string name_;
            Between<T> between_;
    };
    class BaseVariable;
    class Product {
        public:
            Product(){}
            Product(std::string name) : name_{name} {}

            Product(std::string name, const std::initializer_list<Attribute>&& attributes) : name_{name}, attributes_{attributes} {}
            Product(std::string name, const std::set<Attribute>& attributes) : name_{name}, attributes_{attributes} {}

            Product(std::string name, const std::set<Attribute>&& attributes, std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, attributes_{attributes}, conditions_{std::move(conditions)} {}
            Product(std::string name, const std::set<Attribute>& attributes, const std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, attributes_{attributes}, conditions_{conditions} {}

            Product(std::string name, const std::vector<std::shared_ptr<AttributeCondition>>& conditions) : name_{name}, conditions_{conditions} {}
            Product(std::string name, std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, conditions_{std::move(conditions)} {}


            // void name(const std::string& name) {name_.assign(name);}
            void name(std::string name) {name_ = name;}
            const std::string& name() const {return name_;}

            auto& attributes() {return attributes_;}
            auto& conditions() {return conditions_;}

            auto& attributes() const {return attributes_;}
            auto& conditions() const {return conditions_;}

            // void add_attribute(std::unique_ptr<Attribute> attr){attributes_.insert(attributes_.end(), std::move(attr));}
            void add_attribute(const Attribute& attr){attributes_.insert(attributes_.end(), attr);}

            friend bool operator==(const Product& me, const Product& other){
                // return me.name() == other.name() &&
                return boost::iequals(me.name(), other.name()) &&
                    std::is_permutation(other.attributes().cbegin(), other.attributes().cend(), me.attributes().cbegin(), me.attributes().cend()) &&
                    std::is_permutation(other.conditions().cbegin(), other.conditions().cend(), me.conditions().cbegin(), me.conditions().cend());
            }
            bool matches(const std::vector<Product>& other) const {
                for (const auto& o : other){
                    if (matches(o)){
                        return true;
                    }
                }
                return false;
            }
            bool matches(const Product& other) const {
                // std::cout << "product.matches\n";
                // std::cout << "other name: " << other.name() << "\n";
                // std::cout << "this name length: " << name().length() << "\n";
                // std::cout << "this name substr: " << name().substr(0, 1) << "\n";
                // std::cout << "this name: " << name() << "\n";
                if (boost::iequals(name(), other.name())){
                // if (name() == other.name()){
                    // std::cout << "names match\n";
                    const auto& match_atts = other.attributes();
                    
                    bool all_matches = std::all_of(match_atts.cbegin(), match_atts.cend(), [this](const auto& o){return attributes_.find(o) != attributes_.end();});
                    // bool all_matches = true;
                    // for (const auto& a : match_atts){
                    //     if (attributes_.find(a) == attributes_.end()){
                    //         std::cout << "Missing " << a << "\n";
                    //         return false;
                    //     }
                    // }
                    if (all_matches){
                        // std::cout << "all attributes matched, checking conditions\n";
                        all_matches = std::all_of(conditions_.cbegin(), conditions_.cend(), [&match_atts](const auto& c){
                            auto i = std::find_if(match_atts.cbegin(), match_atts.cend(), [&c](const auto& o){return c->matches(o);});
                            return i != match_atts.end();
                        });
                    }
                    return all_matches;
                }
                // std::cout << "no match\n";
                return false;
            }
            friend std::ostream& operator<<(std::ostream& out, const Product& product){
                out << product.name_;
                if (product.attributes().size() != 0 || product.conditions().size()){
                    out << "[";
                    bool printed_first = false;
                    for (auto it = product.attributes().cbegin(); it != product.attributes().cend(); ++it){
                        if (printed_first){
                            out << ',';
                        }
                        out << *it;
                        printed_first = true;
                    }
                    for (auto it = product.conditions().cbegin(); it != product.conditions().cend(); ++it){
                        if (printed_first){
                            out << ',';
                        }
                        out << **it;
                        printed_first = true;
                    }
                    out << "]";
                }
                return out;
            }
            void variable(BaseVariable* variable) {variable_ = variable;}
            BaseVariable*  variable() const {return variable_;}

            static size_t skip_passed_char(const std::string& input, size_t pos, char c, bool in_quotes=false){
                while (pos < input.size()){
                    if (input[pos] == c){
                        return pos + 1;
                    } else if (input[pos] == '\\'){
                        pos += 2;
                    } else if (in_quotes){
                        pos++;
                    } else {
                        switch (input[pos]){
                            case '\'':
                                pos = skip_passed_char(input, pos + 1, '\'', true);
                                break;
                            case '"':
                                pos = skip_passed_char(input, pos + 1, '"', true);
                                break;
                            default:
                                pos++;
                        }
                    }
                }
                return std::string::npos;
            }
            static std::vector<Product> parse_list(const std::string& input){
                std::vector<Product> ret{};

                size_t product_start = 0;
                auto current_pos = product_start;
                while (current_pos != std::string::npos && current_pos < input.size()){
                    if (input[current_pos] == ','){
                        ret.push_back(parse(input.substr(product_start, current_pos - product_start)));
                        current_pos++;
                        product_start = current_pos;
                    } else if (input[current_pos] == '\''){
                        current_pos = skip_passed_char(input, current_pos + 1, '\'', true);
                    } else if (input[current_pos] == '"'){
                        current_pos = skip_passed_char(input, current_pos + 1, '"', true);
                    } else if (input[current_pos] == '['){
                        current_pos = skip_passed_char(input, current_pos + 1, ']');
                    } else {
                        current_pos++;
                    }
                }
                if (current_pos != product_start){
                    ret.push_back(parse(input.substr(product_start, current_pos - product_start)));
                }
                return ret;
            }
            static Product parse(const std::string& input){
                auto attribute_bracket = input.find_first_of('[');

                size_t name_start = 0;
                while (input[name_start] == ' '){
                    name_start++;
                }
                size_t name_end = attribute_bracket;
                if (attribute_bracket == std::string::npos){
                    name_end = input.size();
                }
                while (input[name_end - 1] == ' '){
                    name_end--;
                }

                if (attribute_bracket == std::string::npos){
                    if (name_start == 0 && name_end == attribute_bracket){
                        return Product{input};
                    }
                    return Product{input.substr(name_start, name_end - name_start)};
                }
                std::string name(input, name_start, name_end - name_start);

                std::set<Attribute> attributes{};

                auto param_start = attribute_bracket + 1;
                auto current_pos = param_start;
                while (current_pos != std::string::npos && current_pos < input.size()){
                    if (input[current_pos] == ']'){
                        if (current_pos != param_start){
                            attributes.insert(Attribute::parse(input.substr(param_start, current_pos - param_start)));
                        }
                        break;
                    }
                    if (input[current_pos] == ','){
                        attributes.insert(Attribute::parse(input.substr(param_start, current_pos - param_start)));
                        current_pos++;
                        param_start = current_pos;
                    } else if (input[current_pos] == '\''){
                        current_pos = skip_passed_char(input, current_pos + 1, '\'', true);
                    } else if (input[current_pos] == '"'){
                        current_pos = skip_passed_char(input, current_pos + 1, '"', true);
                    } else {
                        current_pos++;
                    }
                }
                return Product{name, attributes};
            }
        private:
            std::string name_{"unspecified"};
            std::set<Attribute> attributes_{};
            std::vector<std::shared_ptr<AttributeCondition>> conditions_{}; // pointer for hierarchy, shared for copy-able
            std::set<Attribute> configuration_{};
            BaseVariable* variable_{nullptr};


            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}
            // static bool attributes_equal( const Attribute& left, const Attribute& right){return left == right;}
            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}
            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}


    };
} // namespace focs

#endif // FOCS_PRODUCT

