import pathlib
import re
import sys

output = pathlib.Path(sys.argv[1])
classes = {}
selected = set()

for index, directory in enumerate(sys.argv[2:]):
    for header in sorted(pathlib.Path(directory).rglob("*.h")):
        text = header.read_text(encoding="utf-8")
        text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
        for name, base in re.findall(
                r"\bclass\s+(?:OPENMM_EXPORT\w*\s+)?(\w+)"
                r"\s*:\s*public\s+(\w+)", text):
            classes[name] = (base, header)
            if index == 0:
                selected.add(name)

def isForceImpl(name):
    if name == "ForceImpl":
        return True
    return name in classes and isForceImpl(classes[name][0])

names = sorted(name for name in selected
               if isForceImpl(name) and name != "RPMDUpdater")
if not names:
    raise RuntimeError("No ForceImpl subclasses found")

lines = ["#include <type_traits>"]
for name in names:
    lines.append('#include "{}"'.format(
        classes[name][1].resolve().as_posix()))

lines.extend(["", "using namespace OpenMM;", ""])
for name in names:
    lines.append(
        'static_assert(std::is_same<'
        'decltype(&{0}::updateContextState), '
        'void ({0}::*)(ContextImpl&, bool&)>::value, '
        '"{0} must declare the two-argument updateContextState()");'
        .format(name))

lines.extend(["", "int main() {", "    return 0;", "}", ""])
text = "\n".join(lines)
if not output.exists() or output.read_text(encoding="utf-8") != text:
    output.write_text(text, encoding="utf-8")
