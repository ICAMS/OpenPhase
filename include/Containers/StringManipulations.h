/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2024
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef STRINGMANIPULATIONS_H
#define STRINGMANIPULATIONS_H

#include <string>
#include <algorithm>

inline void RemoveAllWhitespaces(std::string& InString)                         ///< Removes all whitespaces (string is modified)
{
    InString.erase(std::remove_if( InString.begin(), InString.end(),
     [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}), InString.end() );
}

inline void RemoveLeadingTrailingWhitespaces(std::string& InString)             ///< Removes leading and trailing whitespaces (string is modified)
{
    std::string whitespaces (" \t");
    size_t first = InString.find_first_not_of(whitespaces);
    size_t last = InString.find_last_not_of(whitespaces);
    if (first == last)
    {
        InString = "";
    }
    else
    {
        InString = InString.substr(first, (last-first+1));
    }
}

inline void ReplacePunctuationMarksWithComma(std::string& InString)             ///< Replaces all punctuation marks with comma (string is modified)
{
    InString.erase (std::remove (InString.begin(), InString.end(), ' '), InString.end());
    std::replace(InString.begin(), InString.end(), '/', ',');
    std::replace(InString.begin(), InString.end(), '-', ',');
    std::replace(InString.begin(), InString.end(), '@', ',');
    std::replace(InString.begin(), InString.end(), ':', ',');
    std::replace(InString.begin(), InString.end(), ';', ',');
    std::replace(InString.begin(), InString.end(), '!', ',');
    std::replace(InString.begin(), InString.end(), '?', ',');
    std::replace(InString.begin(), InString.end(), '|', ',');
}

inline bool StringIsEmpty(std::string& InString)                                ///< Returns true if string is empty or contains only whitespaces
{
    std::string whitespaces (" \t");
    size_t position = InString.find_first_not_of(whitespaces);
    bool result = false;
    if(position == std::string::npos)
    {
        result = true;
    }
    return result;
}

inline bool StringIsComment(std::string& InString, std::string comment_symbol = "//")  ///< Returns true if string is a coment starting with a special comment symbol
{
    size_t position_comment = InString.find(comment_symbol);
    size_t position_text    = InString.find_first_not_of(" \t/");
    bool result = false;
    if(position_comment != std::string::npos and
       position_comment < position_text)
    {
        result = true;
    }
    return result;
}

inline bool StringIsInputParameter(const std::string& InString)                 ///< Returns true if string has a format of input parameter starting with "$Key" and containing ":" delimiter for parameter reading
{
    std::string key_symbol ("$");
    std::string delimiter  (":");

    size_t position_key_begin = InString.find(key_symbol);
    size_t position_key_end   = InString.find_first_of(" \t:", position_key_begin);
    size_t position_delimiter = InString.find(delimiter, position_key_begin);
    bool result = false;
    if(position_key_begin != std::string::npos and
       position_key_end   >  position_key_begin + 2 and
       position_delimiter != std::string::npos and
       position_delimiter >= position_key_end )
    {
        result = true;
    }
    return result;
}

inline bool StringIsSectionName(const std::string& InString)                    ///< Returns true if string is a section name with format @Name, where Name has to contain at least one character
 {
    std::string key_symbol ("@");

    size_t position_key_begin = InString.find(key_symbol);
    size_t position_key_end   = InString.find_first_of(" \t\n:", position_key_begin);
    bool result = false;
    if(position_key_begin != std::string::npos and
       position_key_end   >  position_key_begin + 2)
    {
        result = true;
    }
    return result;
}

inline bool StringIsInBrackets(const std::string& InString, std::string begin_bracket,
                                                     std::string end_bracket)   ///< Returns true if string is inside of a specified brackets pair
{
    size_t position_begin_bracket  = InString.find(begin_bracket);
    size_t position_end_bracket    = InString.find(end_bracket, position_begin_bracket + begin_bracket.size());

    bool result = false;
    if(position_begin_bracket != std::string::npos and
       position_end_bracket   != std::string::npos)
    {
        result = true;
    }
    return result;
}

inline void RemoveBrackets(std::string& InString, std::string begin_bracket,
                                           std::string end_bracket)             ///< Removes brackets from string (string is modified).
{
    size_t position_begin_bracket  = InString.find(begin_bracket);
    size_t position_end_bracket    = InString.find(end_bracket, position_begin_bracket + begin_bracket.size());

    if(position_begin_bracket != std::string::npos and
       position_end_bracket   != std::string::npos)
    {
        InString.erase(position_end_bracket, end_bracket.size());
        InString.erase(position_begin_bracket, begin_bracket.size());
    }
}

inline std::string ExtractStringBetweenBrackets(const std::string& InString,
                                               std::string begin_bracket,
                                               std::string end_bracket)         ///< Extracts string within brackets
{
    size_t position_begin_bracket = InString.find(begin_bracket);
    size_t position_end_bracket   = InString.find(end_bracket, position_begin_bracket + begin_bracket.size());

    std::string result;
    if(position_begin_bracket != std::string::npos and
       position_end_bracket   != std::string::npos)
    {
        result = InString.substr(position_begin_bracket + begin_bracket.size(),
          position_end_bracket - position_begin_bracket - begin_bracket.size());
    }

    return result;
}

inline std::vector<std::string> ExtractSubstrings(const std::string& InString,
                                                  std::string delimiter = ",")  ///< Extracts substrings separated by the delimiter
{
    std::vector<std::string> result;
    std::string tmp = InString;
    bool not_at_end = true;
    size_t position = 0;

    RemoveAllWhitespaces(tmp);

    while(not_at_end)
    {
        not_at_end = false;

        size_t next_delimiter = tmp.find(delimiter, position);
        if(next_delimiter != std::string::npos)
        {
            result.push_back(tmp.substr(position, next_delimiter - position));
            position = next_delimiter + 1;
            not_at_end = true;
        }
        else
        {
            result.push_back(tmp.substr(position));
        }
    }

    return result;
}

#endif
