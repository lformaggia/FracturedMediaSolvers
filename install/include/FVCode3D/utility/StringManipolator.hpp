/*!
 *  @file StringManipolator.hpp
 *  @brief Functions to manipulate strings.
 */

#ifndef STRINGMANIPOLATOR_HPP_
#define STRINGMANIPOLATOR_HPP_

#include <sstream>
#include <string>
#include <vector>

namespace FVCode3D
{

//! Convert a string to the desired type (using a stringstream)
/*!
 *  @param str string to convert
 */
template <typename T>
T lexical_cast( const std::string& str )
{
    T var;
    std::istringstream iss;
    iss.str(str);
    iss >> var;
    return var;
} // lexical_cast

//! Convert a string to its upper case version
inline std::string toUpper ( const std::string& _string )
{
    std::string out( _string );
    std::transform( std::begin( out ), std::end( out ), std::begin( out ), ::toupper );
    return out;
} // toUpper

//! Splits a string into tokes
/*!
 * @warning empty tokens are skipped
 * @param s input string
 * @param[out] tokens tokens read from the string
 * @param delim character separating tokens
 */
inline void tokenize(const std::string &s, std::vector<std::string> &tokens, char delim = ' ')
{
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delim))
    {
        if (!token.empty())
        {
            tokens.push_back(token);
        }
    } // while
} // tokenize

struct MyDelimiters : std::ctype<char>
{
    MyDelimiters():
        std::ctype<char>(get_table()) {}

    static mask const* get_table()
    {
        static mask rc[table_size];
        rc['('] = std::ctype_base::space;
        rc[')'] = std::ctype_base::space;
        rc[' '] = std::ctype_base::space;
        rc[';'] = std::ctype_base::space;
        rc['\n'] = std::ctype_base::space;
        return &rc[0];
    }
}; // MyDelimiters

} // namespace FVCode3D

#endif // STRINGMANIPOLATOR_HPP_
