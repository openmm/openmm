#ifndef LEPTON_OPERATION_H_
#define LEPTON_OPERATION_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "windowsIncludes.h"
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <sstream>

namespace Lepton {

class LEPTON_EXPORT Operation {
public:
    enum Id {CONSTANT, VARIABLE, CUSTOM, ADD, SUBTRACT, MULTIPLY, DIVIDE, POWER, NEGATE, SQRT, EXP, LOG,
             SIN, COS, SEC, CSC, TAN, COT, ASIN, ACOS, ATAN};
    virtual std::string getName() const = 0;
    virtual Id getId() const = 0;
    virtual int getNumArguments() const = 0;
    virtual Operation* clone() const = 0;
    virtual double evaluate(double* args, const std::map<std::string, double>& variables) const = 0;
    class Constant;
    class Variable;
    class Custom;
    class Add;
    class Subtract;
    class Multiply;
    class Divide;
    class Power;
    class Negate;
    class Sqrt;
    class Exp;
    class Log;
    class Sin;
    class Cos;
    class Sec;
    class Csc;
    class Tan;
    class Cot;
    class Asin;
    class Acos;
    class Atan;
};

class Operation::Constant : public Operation {
public:
    Constant(double value) : value(value) {
    }
    std::string getName() const {
        std::stringstream name;
        name << value;
        return name.str();
    }
    Id getId() const {
        return CONSTANT;
    }
    int getNumArguments() const {
        return 0;
    }
    Operation* clone() const {
        return new Constant(value);
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return value;
    }
private:
    double value;
};

class Operation::Variable : public Operation {
public:
    Variable(const std::string& name) : name(name) {
    }
    std::string getName() const {
        return name;
    }
    Id getId() const {
        return VARIABLE;
    }
    int getNumArguments() const {
        return 0;
    }
    Operation* clone() const {
        return new Variable(name);
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        std::map<std::string, double>::const_iterator iter = variables.find(name);
        if (iter == variables.end())
            throw std::exception();
        return iter->second;
    }
private:
    std::string name;
};

class Operation::Custom : public Operation {
public:
    Custom(const std::string& name, int arguments) : name(name), arguments(arguments) {
    }
    std::string getName() const {
        return name;
    }
    Id getId() const {
        return CUSTOM;
    }
    int getNumArguments() const {
        return arguments;
    }
    Operation* clone() const {
        return new Custom(name, arguments);
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return 0.0;
    }
private:
    std::string name;
    int arguments;
};

class Operation::Add : public Operation {
public:
    Add() {
    }
    std::string getName() const {
        return "+";
    }
    Id getId() const {
        return ADD;
    }
    int getNumArguments() const {
        return 2;
    }
    Operation* clone() const {
        return new Add();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return args[0]+args[1];
    }
};

class Operation::Subtract : public Operation {
public:
    Subtract() {
    }
    std::string getName() const {
        return "-";
    }
    Id getId() const {
        return SUBTRACT;
    }
    int getNumArguments() const {
        return 2;
    }
    Operation* clone() const {
        return new Subtract();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return args[0]-args[1];
    }
};

class Operation::Multiply : public Operation {
public:
    Multiply() {
    }
    std::string getName() const {
        return "*";
    }
    Id getId() const {
        return MULTIPLY;
    }
    int getNumArguments() const {
        return 2;
    }
    Operation* clone() const {
        return new Multiply();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return args[0]*args[1];
    }
};

class Operation::Divide : public Operation {
public:
    Divide() {
    }
    std::string getName() const {
        return "/";
    }
    Id getId() const {
        return DIVIDE;
    }
    int getNumArguments() const {
        return 2;
    }
    Operation* clone() const {
        return new Divide();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return args[0]/args[1];
    }
};

class Operation::Power : public Operation {
public:
    Power() {
    }
    std::string getName() const {
        return "^";
    }
    Id getId() const {
        return POWER;
    }
    int getNumArguments() const {
        return 2;
    }
    Operation* clone() const {
        return new Power();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::pow(args[0], args[1]);
    }
};

class Operation::Negate : public Operation {
public:
    Negate() {
    }
    std::string getName() const {
        return "-";
    }
    Id getId() const {
        return NEGATE;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Negate();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return -args[0];
    }
};

class Operation::Sqrt : public Operation {
public:
    Sqrt() {
    }
    std::string getName() const {
        return "sqrt";
    }
    Id getId() const {
        return SQRT;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Sqrt();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::sqrt(args[0]);
    }
};

class Operation::Exp : public Operation {
public:
    Exp() {
    }
    std::string getName() const {
        return "exp";
    }
    Id getId() const {
        return EXP;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Exp();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::exp(args[0]);
    }
};

class Operation::Log : public Operation {
public:
    Log() {
    }
    std::string getName() const {
        return "log";
    }
    Id getId() const {
        return SQRT;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Log();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::log(args[0]);
    }
};

class Operation::Sin : public Operation {
public:
    Sin() {
    }
    std::string getName() const {
        return "sin";
    }
    Id getId() const {
        return LOG;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Sin();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::sin(args[0]);
    }
};

class Operation::Cos : public Operation {
public:
    Cos() {
    }
    std::string getName() const {
        return "cos";
    }
    Id getId() const {
        return COS;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Cos();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::cos(args[0]);
    }
};

class Operation::Sec : public Operation {
public:
    Sec() {
    }
    std::string getName() const {
        return "sec";
    }
    Id getId() const {
        return SEC;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Sec();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return 1.0/std::cos(args[0]);
    }
};

class Operation::Csc : public Operation {
public:
    Csc() {
    }
    std::string getName() const {
        return "csc";
    }
    Id getId() const {
        return CSC;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Csc();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return 1.0/std::sin(args[0]);
    }
};

class Operation::Tan : public Operation {
public:
    Tan() {
    }
    std::string getName() const {
        return "tan";
    }
    Id getId() const {
        return TAN;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Tan();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::tan(args[0]);
    }
};

class Operation::Cot : public Operation {
public:
    Cot() {
    }
    std::string getName() const {
        return "cot";
    }
    Id getId() const {
        return COT;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Cot();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return 1.0/std::tan(args[0]);
    }
};

class Operation::Asin : public Operation {
public:
    Asin() {
    }
    std::string getName() const {
        return "asin";
    }
    Id getId() const {
        return ASIN;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Asin();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::asin(args[0]);
    }
};

class Operation::Acos : public Operation {
public:
    Acos() {
    }
    std::string getName() const {
        return "acos";
    }
    Id getId() const {
        return ACOS;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Acos();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::acos(args[0]);
    }
};

class Operation::Atan : public Operation {
public:
    Atan() {
    }
    std::string getName() const {
        return "atan";
    }
    Id getId() const {
        return ATAN;
    }
    int getNumArguments() const {
        return 1;
    }
    Operation* clone() const {
        return new Atan();
    }
    double evaluate(double* args, const std::map<std::string, double>& variables) const {
        return std::atan(args[0]);
    }
};

} // namespace Lepton

#endif /*LEPTON_OPERATION_H_*/
