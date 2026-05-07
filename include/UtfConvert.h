/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2025
 *  Main contributors :   Raphael Schiedung
 *
 */

#ifndef UTFCONVERT_H
#define UTFCONVERT_H
#include <string>
#include <stdexcept>
#include <cstdint>

#ifdef _WIN32
    #include <windows.h>
#endif

namespace openphase
{

// ================================================================
// Utility: Throw descriptive UTF errors
// ================================================================
inline void throw_unicode_error(const char* msg) {
    throw std::runtime_error(std::string("Unicode conversion error: ") + msg);
}

// ================================================================
// Helper: Determine UTF-8 sequence length from first byte
// ================================================================
inline size_t utf8_sequence_length(uint8_t lead)
{
    if (lead <= 0x7F) return 1;           // ASCII
    if ((lead & 0xE0) == 0xC0) return 2;  // 110xxxxx
    if ((lead & 0xF0) == 0xE0) return 3;  // 1110xxxx
    if ((lead & 0xF8) == 0xF0) return 4;  // 11110xxx
    return 0;                             // Invalid lead byte
}

// ================================================================
// UTF-8 → Unicode codepoint (safe decode)
// Returns: decoded codepoint, advances 'i' by number of bytes
// ================================================================
inline uint32_t decode_utf8_codepoint(const std::string& s, size_t& i)
{
    uint8_t lead = (uint8_t)s[i];
    size_t len = utf8_sequence_length(lead);

    if (len == 0)
        throw_unicode_error("Invalid UTF-8 leading byte");

    if (i + len > s.size())
        throw_unicode_error("Truncated UTF-8 sequence");

    // Decode first byte
    uint32_t cp = lead & ((1 << (8 - len - 1)) - 1);

    // Decode continuation bytes
    for (size_t j = 1; j < len; ++j)
    {
        uint8_t ch = (uint8_t)s[i + j];
        if ((ch & 0xC0) != 0x80)
            throw_unicode_error("Invalid UTF-8 continuation byte");

        cp = (cp << 6) | (ch & 0x3F);
    }

    // Reject overlong sequences (e.g., "\xC0\xAF" for "/")
    static const uint32_t min_values[5] = {
        0, 0, 0x80, 0x800, 0x10000
    };
    if (cp < min_values[len])
        throw_unicode_error("Overlong UTF-8 sequence");

    // Reject surrogate codepoints (U+D800–DFFF)
    if (cp >= 0xD800 && cp <= 0xDFFF)
        throw_unicode_error("UTF-8 sequence encodes invalid surrogate");

    // Reject out-of-range (> U+10FFFF)
    if (cp > 0x10FFFF)
        throw_unicode_error("Invalid Unicode codepoint (> U+10FFFF)");

    i += len;
    return cp;
}

// ================================================================
// Unicode codepoint → UTF-8 (safe encode)
// ================================================================
inline void encode_utf8_codepoint(uint32_t cp, std::string& out)
{
    if (cp <= 0x7F) {
        out.push_back((char)cp);
    }
    else if (cp <= 0x7FF) {
        out.push_back((char)(0xC0 | (cp >> 6)));
        out.push_back((char)(0x80 | (cp & 0x3F)));
    }
    else if (cp <= 0xFFFF) {
        if (cp >= 0xD800 && cp <= 0xDFFF)
            throw_unicode_error("Attempt to encode surrogate half directly");

        out.push_back((char)(0xE0 | (cp >> 12)));
        out.push_back((char)(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back((char)(0x80 | (cp & 0x3F)));
    }
    else if (cp <= 0x10FFFF) {
        out.push_back((char)(0xF0 | (cp >> 18)));
        out.push_back((char)(0x80 | ((cp >> 12) & 0x3F)));
        out.push_back((char)(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back((char)(0x80 | (cp & 0x3F)));
    }
    else {
        throw_unicode_error("Unicode codepoint out of range");
    }
}

// ================================================================
// UTF-8 → UTF-16 (wstring)
// ================================================================
inline std::wstring utf8_to_utf16(const std::string& utf8)
{
#ifdef _WIN32
    // Windows uses UTF-16 wchar_t (16-bit)
    if (utf8.empty()) return {};

    int size_needed = MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS,
        utf8.data(), (int)utf8.size(),
        nullptr, 0
    );

    if (size_needed == 0)
        throw_unicode_error("Invalid UTF-8 input (WinAPI stage 1)");

    std::wstring result(size_needed, 0);

    int converted = MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS,
        utf8.data(), (int)utf8.size(),
        &result[0], size_needed
    );

    if (converted == 0)
        throw_unicode_error("Invalid UTF-8 input (WinAPI stage 2)");

    return result;

#else
    // POSIX: wchar_t is UTF-32 (32-bit), so it's simpler
    std::wstring out;
    out.reserve(utf8.size());

    size_t i = 0;
    while (i < utf8.size()) {
        uint32_t cp = decode_utf8_codepoint(utf8, i);
        out.push_back((wchar_t)cp);
    }

    return out;
#endif
}

// ================================================================
// UTF-16 (wstring) → UTF-8
// ================================================================
inline std::string utf16_to_utf8(const std::wstring& wide)
{
#ifdef _WIN32
    if (wide.empty()) return {};

    int size_needed = WideCharToMultiByte(
        CP_UTF8, WC_ERR_INVALID_CHARS,
        wide.data(), (int)wide.size(),
        nullptr, 0, nullptr, nullptr
    );

    if (size_needed == 0)
        throw_unicode_error("Invalid UTF-16 input (WinAPI stage 1)");

    std::string result(size_needed, 0);

    int converted = WideCharToMultiByte(
        CP_UTF8, WC_ERR_INVALID_CHARS,
        wide.data(), (int)wide.size(),
        &result[0], size_needed, nullptr, nullptr
    );

    if (converted == 0)
        throw_unicode_error("Invalid UTF-16 input (WinAPI stage 2)");

    return result;

#else
    // Linux/macOS: wchar_t is 32-bit, holding full codepoints
    std::string out;
    out.reserve(wide.size());

    for (wchar_t wc : wide) {
        uint32_t cp = (uint32_t)wc;
        encode_utf8_codepoint(cp, out);
    }

    return out;
#endif
}
#endif // UTFCONVERT_H
} // namespace openphase
