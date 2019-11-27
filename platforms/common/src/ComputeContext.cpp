/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/common/ComputeContext.h"
#include <set>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

string ComputeContext::replaceStrings(const string& input, const std::map<std::string, std::string>& replacements) const {
    static set<char> symbolChars;
    if (symbolChars.size() == 0) {
        symbolChars.insert('_');
        for (char c = 'a'; c <= 'z'; c++)
            symbolChars.insert(c);
        for (char c = 'A'; c <= 'Z'; c++)
            symbolChars.insert(c);
        for (char c = '0'; c <= '9'; c++)
            symbolChars.insert(c);
    }
    string result = input;
    for (auto& pair : replacements) {
        int index = 0;
        int size = pair.first.size();
        do {
            index = result.find(pair.first, index);
            if (index != result.npos) {
                if ((index == 0 || symbolChars.find(result[index-1]) == symbolChars.end()) && (index == result.size()-size || symbolChars.find(result[index+size]) == symbolChars.end())) {
                    // We have found a complete symbol, not part of a longer symbol.

                    result.replace(index, size, pair.second);
                    index += pair.second.size();
                }
                else
                    index++;
            }
        } while (index != result.npos);
    }
    return result;
}

string ComputeContext::doubleToString(double value) const {
    stringstream s;
    s.precision(getUseDoublePrecision() ? 16 : 8);
    s << scientific << value;
    if (!getUseDoublePrecision())
        s << "f";
    return s.str();
}

string ComputeContext::intToString(int value) const {
    stringstream s;
    s << value;
    return s.str();
}
