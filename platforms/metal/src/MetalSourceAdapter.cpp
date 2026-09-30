/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalSourceAdapter.h"
#include "MetalPairwiseOptimizations.h"
#include "MetalReductionOptimizations.h"
#include "CommonKernelSources.h"
#include "openmm/OpenMMException.h"
#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace std;

namespace {

struct Token {
    string text;
    size_t begin, end;
    bool directive;
};

bool identifier(char c) {
    return isalnum(static_cast<unsigned char>(c)) || c == '_';
}

/** Lex only the structure needed by the adapter; retain source offsets/comments. */
vector<Token> tokenize(const string& source) {
    vector<Token> tokens;
    for (size_t i = 0; i < source.size();) {
        if (isspace(static_cast<unsigned char>(source[i]))) { i++; continue; }
        if (source.compare(i, 2, "//") == 0) {
            size_t end = source.find('\n', i);
            i = (end == string::npos ? source.size() : end);
            continue;
        }
        if (source.compare(i, 2, "/*") == 0) {
            size_t end = source.find("*/", i+2);
            if (end == string::npos)
                throw OpenMMException("Unterminated comment in Metal Common source");
            i = end+2;
            continue;
        }
        size_t begin = i;
        bool directive = source[i] == '#';
        if (directive) {
            do {
                size_t end = source.find('\n', i);
                i = (end == string::npos ? source.size() : end+1);
            } while (i < source.size() && i >= 2 && source[i-2] == '\\');
        }
        else if (source[i] == '"' || source[i] == '\'') {
            char quote = source[i++];
            while (i < source.size()) {
                if (source[i++] == '\\') { if (i < source.size()) i++; }
                else if (source[i-1] == quote) break;
            }
        }
        else if (identifier(source[i])) {
            while (i < source.size() && identifier(source[i])) i++;
        }
        else i++;
        tokens.push_back({source.substr(begin, i-begin), begin, i, directive});
    }
    return tokens;
}

size_t matching(const vector<Token>& tokens, size_t start, const string& left, const string& right) {
    int depth = 0;
    struct Conditional { int initialDepth, branchDepth; bool hasAlternative; };
    vector<Conditional> conditionals;
    for (size_t i = start; i < tokens.size(); i++) {
        if (tokens[i].directive) {
            istringstream directive(tokens[i].text.substr(1));
            string keyword;
            directive >> keyword;
            if (keyword == "if" || keyword == "ifdef" || keyword == "ifndef")
                conditionals.push_back({depth, depth, false});
            else if (!conditionals.empty() && (keyword == "else" || keyword == "elif")) {
                Conditional& branch = conditionals.back();
                if (branch.hasAlternative && depth != branch.branchDepth)
                    throw OpenMMException("Inconsistent conditional delimiters in Metal Common source");
                branch.branchDepth = depth;
                branch.hasAlternative = true;
                depth = branch.initialDepth;
            }
            else if (!conditionals.empty() && keyword == "endif") {
                const Conditional& branch = conditionals.back();
                if (branch.hasAlternative && depth != branch.branchDepth)
                    throw OpenMMException("Inconsistent conditional delimiters in Metal Common source");
                conditionals.pop_back();
            }
            continue;
        }
        if (tokens[i].text == left) depth++;
        if (tokens[i].text == right && --depth == 0) return i;
    }
    throw OpenMMException("Unbalanced delimiters in Metal Common source");
}

struct Function {
    string name;
    size_t start, nameToken, open, close, body, end;
    bool kernel;
};

vector<Function> findFunctions(const vector<Token>& tokens) {
    vector<Function> functions;
    size_t statement = 0;
    for (size_t i = 0; i < tokens.size(); i++) {
        if (tokens[i].directive) { statement = i+1; continue; }
        if (tokens[i].text == ";" || tokens[i].text == "}") { statement = i+1; continue; }
        if (tokens[i].text == "{" ) {
            i = matching(tokens, i, "{", "}");
            statement = i+1;
            continue;
        }
        if (tokens[i].text != "(" || i == 0 || !identifier(tokens[i-1].text[0])) continue;
        size_t close = matching(tokens, i, "(", ")");
        if (close+1 >= tokens.size() || tokens[close+1].text != "{") { i = close; continue; }
        bool kernel = false;
        for (size_t j = statement; j < i; j++)
            kernel |= (tokens[j].text == "KERNEL" || tokens[j].text == "__kernel");
        size_t end = matching(tokens, close+1, "{", "}");
        functions.push_back({tokens[i-1].text, statement, i-1, i, close, close+1, end, kernel});
        i = end;
        statement = i+1;
    }
    return functions;
}

struct Argument {
    string declaration, name, directive;
    bool local;
};

vector<Argument> parseArguments(const vector<Token>& tokens, size_t begin, size_t end) {
    vector<Argument> result;
    vector<string> words;
    auto finish = [&]() {
        if (words.empty()) return;
        if (words.size() == 1 && words[0] == "void") { words.clear(); return; }
        string declaration, name = words.back();
        bool local = false;
        if (name.empty() || !identifier(name[0]) || name == "*")
            throw OpenMMException("Unsupported Metal Common kernel argument declaration: "+name);
        for (const string& word : words) {
            if (word == "volatile" || word == "RESTRICT" || word == "restrict") continue;
            local |= (word == "LOCAL_ARG" || word == "__local");
            declaration += word+" ";
        }
        result.push_back({declaration, name, "", local});
        words.clear();
    };
    for (size_t i = begin; i < end; i++) {
        if (tokens[i].directive) {
            finish();
            result.push_back({"", "", tokens[i].text, false});
        }
        else if (tokens[i].text == ",") finish();
        else words.push_back(tokens[i].text);
    }
    finish();
    return result;
}

struct Edit { size_t begin, end; string replacement; };

/** Extract a complete, unchanged function from the corresponding Common template. */
string templateFunction(const string& source, const string& name) {
    const vector<Token> tokens = tokenize(source);
    for (const Function& function : findFunctions(tokens))
        if (function.name == name)
            return source.substr(tokens[function.start].begin, tokens[function.end].end-tokens[function.start].begin);
    throw OpenMMException("Missing Common function for a Metal fast path: "+name);
}

/**
 * Select two existing CUDA branches only inside exact Common template functions.
 * Keep mathematical bodies and synchronization unchanged; never impersonate CUDA
 * or HIP for the rest of a program. The ballot kernel visits padded 32-lane warps.
 */
void appendFastPathEdits(const string& source, const vector<Token>& tokens,
        const vector<Function>& functions, vector<Edit>& edits) {
    for (const Function& function : functions) {
        if (function.name != "reduceMax" && function.name != "findNeighbors") continue;
        bool shuffle = false, ballot = false;
        const string body = source.substr(tokens[function.start].begin,
                tokens[function.end].end-tokens[function.start].begin);
#if OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE
        if (function.name == "reduceMax") {
            static const string original = templateFunction(CommonKernelSources::customNonbondedGroups, "reduceMax");
            shuffle = (body == original);
        }
#endif
#if OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_BALLOT
        if (function.name == "findNeighbors") {
            static const string original = templateFunction(CommonKernelSources::customManyParticle, "findNeighbors");
            ballot = (body == original);
        }
#endif
        if (!shuffle && !ballot) continue;
        int selectedBranches = 0, skippedBranches = 0, calls = 0, scans = 0, masks = 0;
        for (size_t i = function.body+1; i < function.end; i++) {
            const Token& token = tokens[i];
            if (token.directive) {
                string directive;
                for (char c : token.text)
                    if (!isspace(static_cast<unsigned char>(c))) directive += c;
                if ((shuffle && directive == "#ifdefined(__CUDA_ARCH__)&&__CUDA_ARCH__>=700") ||
                        (ballot && directive == "#ifdefined(__CUDA_ARCH__)||defined(USE_HIP)")) {
                    edits.push_back({token.begin, token.end, "#if 1 // Scoped Metal CUDA-derived fast path\n"});
                    selectedBranches++;
                }
                else if (ballot && directive == "#if!(defined(__CUDA_ARCH__)||defined(USE_HIP))") {
                    edits.push_back({token.begin, token.end, "#if 0 // Ballot replaces the local flag array\n"});
                    skippedBranches++;
                }
                continue;
            }
            if (shuffle && token.text == "__shfl_xor_sync" && i+3 < function.end &&
                    tokens[i+1].text == "(" && tokens[i+2].text == "0xffffffff" && tokens[i+3].text == ",") {
                edits.push_back({token.begin, tokens[i+3].end, "simd_shuffle_xor("});
                calls++;
            }
            if (ballot && (token.text == "BALLOT" || token.text == "__ffs") && tokens[i+1].text == "(") {
                size_t close = matching(tokens, i+1, "(", ")");
                if (token.text == "BALLOT") {
                    edits.push_back({token.begin, token.end, "uint(simd_vote::vote_t(simd_ballot"});
                    edits.push_back({tokens[close].end, tokens[close].end, "))"});
                    calls++;
                }
                else {
                    // The surrounding while guarantees a nonzero mask, so
                    // CUDA's one-based ffs is precisely ctz(mask)+1 here.
                    edits.push_back({token.begin, token.end, "(ctz"});
                    edits.push_back({tokens[close].end, tokens[close].end, "+1)"});
                    scans++;
                }
            }
            if (ballot && token.text == "int" && tokens[i+1].text == "includeBlockFlags") {
                // A lane-31-only ballot is INT_MIN as a signed integer; its
                // flags-1 step must use unsigned modulo arithmetic.
                edits.push_back({token.begin, token.end, "uint"});
                masks++;
            }
        }
        if (selectedBranches != 1 || calls != 1 ||
                (ballot && (skippedBranches != 1 || scans != 1 || masks != 1)))
            throw OpenMMException("Common template changed: review the Metal fast path for "+function.name);
    }
}

/**
 * Rewrite the audited fixed-point ABI for scoped floating-accumulator execution.
 * Integer tile counts, matrix indices, and the DPD random state stay 64-bit.
 * This pass precedes signature adaptation, so reflection sees float pointers.
 */
string floatingAccumulatorSource(const string& source, bool inspectFunctions=true) {
    const vector<Token> tokens = tokenize(source);
    const vector<Function> functions = inspectFunctions ? findFunctions(tokens) : vector<Function>();
    const set<string> forceScalars = {"fx", "fy", "fz", "fx0", "fy0", "fz0", "fx1", "fy1", "fz1", "zero"};
    const set<string> integerCastOperands = {"numTiles", "NUM_BLOCKS", "NUM_TILES", "ii", "jj"};
    vector<Edit> edits;
    for (size_t i = 0; i < tokens.size(); i++) {
        const Token& token = tokens[i];
        if (token.directive) {
            // Generated CustomGB derivative macros contain fixed-point casts.
            edits.push_back({token.begin, token.end, "#"+floatingAccumulatorSource(token.text.substr(1), false)});
            continue;
        }
        string functionName;
        for (const Function& function : functions)
            if (i >= function.start && i <= function.end) { functionName = function.name; break; }
        // An inactive ATM state can have an overlapping particle pair and an
        // infinite force. Its exactly zero weight must contribute zero, not
        // IEEE 0*Inf (NaN), to the active state's floating accumulator.
        if (functionName == "hybridForce" && (token.text == "dEdu0" || token.text == "dEdu1") &&
                i+2 < tokens.size() && tokens[i+1].text == "*" && forceScalars.count(tokens[i+2].text))
            edits.push_back({token.begin, tokens[i+2].end,
                "("+token.text+" == 0 ? 0.0f : "+token.text+"*"+tokens[i+2].text+")"});
        if ((token.text == "0x100000000" && functionName != "getRandomNormal") ||
                (token.text == "0xFFFFFFFF" && functionName == "computePerDof"))
            edits.push_back({token.begin, token.end, "1.0f"});
        if (functionName == "convertForces" && token.text == "if" && i+1 < tokens.size() && tokens[i+1].text == "(") {
            size_t end = matching(tokens, i+1, "(", ")");
            bool checksLimit = false;
            for (size_t j = i+2; j < end; j++) checksLimit |= tokens[j].text == "limit";
            if (checksLimit)
                edits.push_back({tokens[i+1].end, tokens[end].begin, "!isfinite(fx) || !isfinite(fy) || !isfinite(fz)"});
        }
        if (token.text != "mm_long" && token.text != "mm_ulong" && token.text != "long" && token.text != "ulong")
            continue;
        bool convert = i+1 < tokens.size() && tokens[i+1].text == "*";
        if (i+1 < tokens.size()) {
            const string& name = tokens[i+1].text;
            // GBSA's generated AtomData member and per-particle locals also
            // hold fixed-point Born derivatives, despite not being pointers.
            convert |= forceScalars.count(name) || name.find("bornForce") != string::npos;
        }
        if (i > 0 && i+1 < tokens.size() && tokens[i-1].text == "(" && tokens[i+1].text == ")") {
            size_t operand = i+2;
            while (operand < tokens.size() && (tokens[operand].text == "(" || tokens[operand].text == "-" || tokens[operand].text == "+")) operand++;
            if (operand < tokens.size() && !integerCastOperands.count(tokens[operand].text)) {
                const string& value = tokens[operand].text;
                convert = forceScalars.count(value) || value.compare(0, 5, "force") == 0 ||
                    value == "realToFixedPoint" || value == "kdr" || value == "dEdu0" ||
                    value == "sum" || value == "mm_long" || value == "mm_ulong";
                if (!convert)
                    throw OpenMMException("Unaudited 64-bit cast in Metal floating accumulators: "+value);
            }
        }
        if (convert) {
            size_t begin = token.begin;
            if (token.text == "long" && i > 0 && tokens[i-1].text == "unsigned") begin = tokens[i-1].begin;
            edits.push_back({begin, token.end, "float"});
        }
    }
    sort(edits.begin(), edits.end(), [](const Edit& a, const Edit& b) { return a.begin < b.begin; });
    string result;
    size_t previous = 0;
    for (const Edit& edit : edits) {
        if (edit.begin < previous)
            throw OpenMMException("Overlapping Metal floating-accumulator adaptations");
        result += source.substr(previous, edit.begin-previous)+edit.replacement;
        previous = edit.end;
    }
    return result+source.substr(previous);
}

bool vectorType(const string& type) {
    static const set<string> types = {"real2", "real3", "real4", "mixed2", "mixed3", "mixed4",
        "float2", "float3", "float4", "int2", "int3", "int4", "uint2", "uint3", "uint4",
        "short2", "short3", "short4", "ushort2", "ushort3", "ushort4", "long2", "long3", "long4"};
    return types.count(type) != 0;
}

string kernelPrefix(const Function& function, const vector<Argument>& arguments) {
    const string prefix = "_metal_"+function.name;
    string ids = "enum { "+prefix+"_base = __COUNTER__ };\nenum {\n";
    string fields, locals, localArguments;
    for (const auto& arg : arguments) {
        if (!arg.directive.empty()) {
            ids += "\n"+arg.directive+"\n";
            fields += "\n"+arg.directive+"\n";
            locals += "\n"+arg.directive+"\n";
            localArguments += "\n"+arg.directive+"\n";
            continue;
        }
        const string id = prefix+"_id_"+arg.name;
        ids += id+" = __COUNTER__-"+prefix+"_base-1,\n";
        if (arg.local) {
            fields += "uint _metal_local_"+arg.name+" [[id("+id+")]];\n";
            localArguments += ", "+arg.declaration+" [[threadgroup("+id+")]]\n";
        }
        else {
            fields += arg.declaration+" [[id("+id+")]];\n";
            locals += "auto "+arg.name+" = _metal_args."+arg.name+";\n";
        }
    }
    ids += prefix+"_count = __COUNTER__-"+prefix+"_base-1\n};\n";
    if (fields.empty()) fields = "uint _metal_unused [[id(0)]];\n";
    return ids+"struct "+prefix+"_arguments {\n"+fields+"};\n"
        "kernel void "+function.name+"(constant "+prefix+"_arguments& _metal_args [[buffer(0)]],\n"
        "uint _metal_gid [[thread_position_in_grid]], uint _metal_lid [[thread_position_in_threadgroup]],\n"
        "uint _metal_group [[threadgroup_position_in_grid]], uint _metal_size [[threads_per_threadgroup]],\n"
        "uint _metal_groups [[threadgroups_per_grid]]\n"+localArguments+
        "#if OPENMM_METAL_CHECK_FIXED_POINT_RANGE\n"
        ", device atomic_uint* _metal_fixedPointRange [[buffer(1)]]\n"
        "#endif\n) {\n"
        "const MetalExecutionContext _metal = {_metal_gid, _metal_lid, _metal_group, _metal_size, _metal_groups\n"
        "#if OPENMM_METAL_CHECK_FIXED_POINT_RANGE\n, _metal_fixedPointRange\n#endif\n};\n"+locals;
}

} // namespace

bool MetalSourceAdapter::isCommonSource(const string& source) {
    for (const Token& token : tokenize(source))
        if (!token.directive && (token.text == "KERNEL" || token.text == "__kernel")) return true;
    return false;
}

string MetalSourceAdapter::translate(const string& originalSource, bool floatingAccumulators) {
    // Optimize the unmodified Common templates once, before changing their
    // accumulator representation or adding Metal address spaces and bindings.
    string source = MetalPairwiseOptimizations::apply(MetalReductionOptimizations::apply(originalSource));
    if (floatingAccumulators)
        source = floatingAccumulatorSource(source);
    const vector<Token> tokens = tokenize(source);
    const vector<Function> functions = findFunctions(tokens);
    vector<Edit> edits;
    appendFastPathEdits(source, tokens, functions, edits);
    set<string> helpers;
    set<size_t> declarations;
    vector<pair<size_t, size_t>> replaced;
    for (const auto& function : functions) {
        declarations.insert(function.nameToken);
        if (function.kernel) {
            edits.push_back({tokens[function.start].begin, tokens[function.body].end,
                kernelPrefix(function, parseArguments(tokens, function.open+1, function.close))});
            replaced.push_back({function.start, function.body});
        }
        else {
            helpers.insert(function.name);
            const string separator = (function.close == function.open+1 ? "" : ", ");
            edits.push_back({tokens[function.open].begin, tokens[function.open].end,
                "(MetalExecutionContext _metal"+separator});
            // OpenCL private helper pointers become explicitly thread-qualified.
            size_t start = function.open+1;
            int nesting = 0;
            for (size_t i = start; i <= function.close; i++) {
                if (tokens[i].text == "(") nesting++;
                if (tokens[i].text == ")" && i != function.close) nesting--;
                if (i != function.close && (tokens[i].text != "," || nesting != 0)) continue;
                bool pointer = false, qualified = false;
                size_t first = start;
                while (first < i && tokens[first].directive) first++;
                for (size_t j = first; j < i; j++) {
                    pointer |= tokens[j].text == "*" || tokens[j].text == "[";
                    qualified |= tokens[j].text == "GLOBAL" || tokens[j].text == "LOCAL_ARG" ||
                        tokens[j].text == "__global" || tokens[j].text == "__local" ||
                        tokens[j].text == "thread" || tokens[j].text == "device" || tokens[j].text == "constant";
                }
                if (pointer && !qualified && first < i)
                    edits.push_back({tokens[first].begin, tokens[first].begin, "thread "});
                start = i+1;
            }
            // Private pointer-to-array locals occur in the Gay-Berne matrix
            // helpers. Their pointee remains in the calling thread's storage.
            for (size_t i = function.body+1; i+5 < function.end; i++) {
                if (identifier(tokens[i].text[0]) && tokens[i+1].text == "(" && tokens[i+2].text == "*" &&
                        identifier(tokens[i+3].text[0]) && tokens[i+4].text == ")" && tokens[i+5].text == "[" &&
                        tokens[i-1].text != "GLOBAL" && tokens[i-1].text != "LOCAL_ARG" &&
                        tokens[i-1].text != "__global" && tokens[i-1].text != "__local")
                    edits.push_back({tokens[i].begin, tokens[i].begin, "thread "});
            }
            // Route the Common helper through the selected Metal primitive:
            // OpenCL-style CAS with an atomic initial load by default, or
            // native float addition under its independent build switch.
            if (function.name == "atomicAddMixed") {
                static const string original = templateFunction(CommonKernelSources::minimize, "atomicAddMixed");
                static const string floating = floatingAccumulatorSource(original);
                const string helper = source.substr(tokens[function.start].begin,
                        tokens[function.end].end-tokens[function.start].begin);
                if (helper == original || helper == floating) {
                    edits.push_back({tokens[function.body].end, tokens[function.end].begin,
                        "\nmetalAtomicAdd(target, value);\n"});
                    replaced.push_back({function.body+1, function.end-1});
                }
            }
        }
    }
    for (size_t i = 0; i < tokens.size(); i++) {
        if (tokens[i].directive) {
            // Macro bodies can call the same DEVICE helpers as ordinary code.
            // Preserve the directive itself and only thread through those calls.
            vector<Token> directive = tokenize(tokens[i].text.substr(1));
            if (!directive.empty() && directive[0].text == "define") {
                for (size_t j = 2; j+1 < directive.size(); j++) {
                    if (helpers.count(directive[j].text) && directive[j+1].text == "(") {
                        size_t position = tokens[i].begin+1+directive[j+1].end;
                        edits.push_back({position, position,
                            (j+2 < directive.size() && directive[j+2].text == ")" ? "_metal" : "_metal, ")});
                    }
                }
            }
        }
        bool skip = tokens[i].directive;
        for (const auto& range : replaced) skip |= i >= range.first && i <= range.second;
        if (skip) continue;
        if (tokens[i].text == "thread") edits.push_back({tokens[i].begin, tokens[i].end, "_metal_thread"});
        if (tokens[i].text == "volatile") edits.push_back({tokens[i].begin, tokens[i].end, ""});
        if (tokens[i].text.size() > 3 && tokens[i].text.substr(tokens[i].text.size()-3) == "ULL")
            edits.push_back({tokens[i].end-3, tokens[i].end, "UL"});
        if (i+1 < tokens.size() && tokens[i+1].text == "(" && helpers.count(tokens[i].text) && !declarations.count(i))
            edits.push_back({tokens[i+1].end, tokens[i+1].end,
                (i+2 < tokens.size() && tokens[i+2].text == ")" ? "_metal" : "_metal, ")});
        if (i+3 < tokens.size() && tokens[i].text == "(" && vectorType(tokens[i+1].text) &&
                tokens[i+2].text == ")" && tokens[i+3].text == "(") {
            edits.push_back({tokens[i].begin, tokens[i+2].end, tokens[i+1].text});
            i += 2;
        }
    }
    sort(edits.begin(), edits.end(), [](const Edit& a, const Edit& b) {
        return a.begin < b.begin || (a.begin == b.begin && a.end < b.end);
    });
    string result;
    size_t previous = 0;
    for (const Edit& edit : edits) {
        if (edit.begin < previous)
            throw OpenMMException("Overlapping Metal Common source adaptations");
        result += source.substr(previous, edit.begin-previous)+edit.replacement;
        previous = edit.end;
    }
    return result+source.substr(previous);
}
