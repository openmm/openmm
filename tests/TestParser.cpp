#include "../libraries/lepton/include/Lepton.h"
#include "openmm/internal/AssertionUtilities.h"

#include <iostream>
#include <limits>
#include <map>

using namespace Lepton;
using namespace OpenMM;
using namespace std;

/**
 * This is a custom function equal to f(x,y) = 2*x*y.
 */

class ExampleFunction : public CustomFunction {
    int getNumArguments() const {
        return 2;
    }
    double evaluate(const double* arguments) const {
        return 2.0*arguments[0]*arguments[1];
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        if (derivOrder[0] == 1) {
            if (derivOrder[1] == 0)
                return 2.0*arguments[1];
            else if (derivOrder[1] == 1)
                return 2.0;
        }
        if (derivOrder[1] == 1 && derivOrder[0] == 0)
            return 2.0*arguments[0];
        return 0.0;
    }
    CustomFunction* clone() const {
        return new ExampleFunction();
    }
};

/**
 * Verify that an expression gives the correct value.
 */

void verifyEvaluation(const string& expression, double expectedValue) {
    map<string, CustomFunction*> customFunctions;
    ParsedExpression parsed = Parser::parse(expression, customFunctions);
    double value = parsed.evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Try optimizing it and make sure the result is still correct.

    value = parsed.optimize().evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Create an ExpressionProgram and see if that also gives the same result.

    ExpressionProgram program = parsed.createProgram();
    value = program.evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Create a CompiledExpression and see if that also gives the same result.

    CompiledExpression compiled = parsed.createCompiledExpression();
    value = compiled.evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);
}

/**
 * Verify that an expression with variables gives the correct value.
 */

void verifyEvaluation(const string& expression, double x, double y, double expectedValue) {
    map<string, double> variables;
    variables["x"] = x;
    variables["y"] = y;
    ParsedExpression parsed = Parser::parse(expression);
    double value = parsed.evaluate(variables);
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Try optimizing it and make sure the result is still correct.

    value = parsed.optimize().evaluate(variables);
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Try optimizing with predefined values for the variables.

    value = parsed.optimize(variables).evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Create an ExpressionProgram and see if that also gives the same result.

    ExpressionProgram program = parsed.createProgram();
    value = program.evaluate(variables);
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);

    // Create a CompiledExpression and see if that also gives the same result.

    CompiledExpression compiled = parsed.createCompiledExpression();
    if (compiled.getVariables().find("x") != compiled.getVariables().end())
        compiled.getVariableReference("x") = x;
    if (compiled.getVariables().find("y") != compiled.getVariables().end())
        compiled.getVariableReference("y") = y;
    value = compiled.evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);
    
    // Try specifying memory locations for the compiled expression.
    
    map<string, double*> variablePointers;
    variablePointers["x"] = &x;
    variablePointers["y"] = &y;
    CompiledExpression compiled2 = parsed.createCompiledExpression();
    compiled2.setVariableLocations(variablePointers);
    value = compiled2.evaluate();
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);
    ASSERT_EQUAL(&x, &compiled2.getVariableReference("x"));
    ASSERT_EQUAL(&y, &compiled2.getVariableReference("y"));

    // Make sure that variable renaming works.

    variables.clear();
    variables["w"] = x;
    variables["y"] = y;
    map<string, string> replacements;
    replacements["x"] = "w";
    value = parsed.renameVariables(replacements).evaluate(variables);
    ASSERT_EQUAL_TOL(expectedValue, value, 1e-10);
}

/**
 * Confirm that a parse error gets thrown.
 */

void verifyInvalidExpression(const string& expression) {
    try {
        Parser::parse(expression);
    }
    catch (const exception& ex) {
        return;
    }
    throw exception();
}

/**
 * Verify that two numbers have the same value.
 */

void assertNumbersEqual(double val1, double val2) {
    const double inf = numeric_limits<double>::infinity();
    if (val1 == val1 || val2 == val2) // If both are NaN, that's fine.
        if (val1 != inf || val2 != inf) // Both infinity is also fine.
            if (val1 != -inf || val2 != -inf) // Same for -infinity.
                ASSERT_EQUAL_TOL(val1, val2, 1e-10);
}

/**
 * Verify that two expressions give the same value.
 */

void verifySameValue(const ParsedExpression& exp1, const ParsedExpression& exp2, double x, double y) {
    map<string, double> variables;
    variables["x"] = x;
    variables["y"] = y;
    double val1 = exp1.evaluate(variables);
    double val2 = exp2.evaluate(variables);
    assertNumbersEqual(val1, val2);
    
    // Now create CompiledExpressions from them and see if those also match.

    CompiledExpression compiled1 = exp1.createCompiledExpression();
    CompiledExpression compiled2 = exp2.createCompiledExpression();
    if (compiled1.getVariables().find("x") != compiled1.getVariables().end())
        compiled1.getVariableReference("x") = x;
    if (compiled1.getVariables().find("y") != compiled1.getVariables().end())
        compiled1.getVariableReference("y") = y;
    if (compiled2.getVariables().find("x") != compiled2.getVariables().end())
        compiled2.getVariableReference("x") = x;
    if (compiled2.getVariables().find("y") != compiled2.getVariables().end())
        compiled2.getVariableReference("y") = y;
    assertNumbersEqual(val1, compiled1.evaluate());
    assertNumbersEqual(val2, compiled2.evaluate());
}

/**
 * Verify that the derivative of an expression is calculated correctly.
 */

void verifyDerivative(const string& expression, const string& expectedDeriv) {
    ParsedExpression computed = Parser::parse(expression).differentiate("x").optimize();
    ParsedExpression expected = Parser::parse(expectedDeriv);
    verifySameValue(computed, expected, 1.0, 2.0);
    verifySameValue(computed, expected, 2.0, 3.0);
    verifySameValue(computed, expected, -2.0, 3.0);
    verifySameValue(computed, expected, 2.0, -3.0);
    verifySameValue(computed, expected, 0.0, -3.0);
    verifySameValue(computed, expected, 2.0, 0.0);
}

/**
 * Test the use of a custom function.
 */

void testCustomFunction(const string& expression, const string& equivalent) {
    map<string, CustomFunction*> functions;
    ExampleFunction exp;
    functions["custom"] = &exp;
    ParsedExpression exp1 = Parser::parse(expression, functions);
    ParsedExpression exp2 = Parser::parse(equivalent);
    verifySameValue(exp1, exp2, 1.0, 2.0);
    verifySameValue(exp1, exp2, 2.0, 3.0);
    verifySameValue(exp1, exp2, -2.0, 3.0);
    verifySameValue(exp1, exp2, 2.0, -3.0);
    ParsedExpression deriv1 = exp1.differentiate("x").optimize();
    ParsedExpression deriv2 = exp2.differentiate("x").optimize();
    verifySameValue(deriv1, deriv2, 1.0, 2.0);
    verifySameValue(deriv1, deriv2, 2.0, 3.0);
    verifySameValue(deriv1, deriv2, -2.0, 3.0);
    verifySameValue(deriv1, deriv2, 2.0, -3.0);
    ParsedExpression deriv3 = deriv1.differentiate("y").optimize();
    ParsedExpression deriv4 = deriv2.differentiate("y").optimize();
    verifySameValue(deriv3, deriv4, 1.0, 2.0);
    verifySameValue(deriv3, deriv4, 2.0, 3.0);
    verifySameValue(deriv3, deriv4, -2.0, 3.0);
    verifySameValue(deriv3, deriv4, 2.0, -3.0);
}

int main() {
    try {
        verifyEvaluation("5", 5.0);
        verifyEvaluation("5*2", 10.0);
        verifyEvaluation("2*3+4*5", 26.0);
        verifyEvaluation("2^-3", 0.125);
        verifyEvaluation("1e+2", 100.0);
        verifyEvaluation("-x", 2.0, 3.0, -2.0);
        verifyEvaluation("y^-x", 3.0, 2.0, 0.125);
        verifyEvaluation("1/-x", 3.0, 2.0, -1.0/3.0);
        verifyEvaluation("2.1e-4*x*(y+1)", 3.0, 1.0, 1.26e-3);
        verifyEvaluation("sin(2.5)", std::sin(2.5));
        verifyEvaluation("cot(x)", 3.0, 1.0, 1.0/std::tan(3.0));
        verifyEvaluation("x^2+y^3+x^-1+y^(1/2)", 1.0, 1.0, 4.0);
        verifyEvaluation("(2*x)*3", 4.0, 4.0, 24.0);
        verifyEvaluation("(x*2)*3", 4.0, 4.0, 24.0);
        verifyEvaluation("2*(x*3)", 4.0, 4.0, 24.0);
        verifyEvaluation("2*(3*x)", 4.0, 4.0, 24.0);
        verifyEvaluation("2*x/3", 1.0, 4.0, 2.0/3.0);
        verifyEvaluation("x*2/3", 1.0, 4.0, 2.0/3.0);
        verifyEvaluation("5*(-x)*(-y)", 1.0, 4.0, 20.0);
        verifyEvaluation("5*(-x)*(y)", 1.0, 4.0, -20.0);
        verifyEvaluation("5*(x)*(-y)", 1.0, 4.0, -20.0);
        verifyEvaluation("5*(-x)/(-y)", 1.0, 4.0, 1.25);
        verifyEvaluation("5*(-x)/(y)", 1.0, 4.0, -1.25);
        verifyEvaluation("5*(x)/(-y)", 1.0, 4.0, -1.25);
        verifyEvaluation("x+(-y)", 1.0, 4.0, -3.0);
        verifyEvaluation("(-x)+y", 1.0, 4.0, 3.0);
        verifyEvaluation("x/(1/y)", 1.0, 4.0, 4.0);
        verifyEvaluation("x*w; w = 5", 3.0, 1.0, 15.0);
        verifyEvaluation("a+b^2;a=x-b;b=3*y", 2.0, 3.0, 74.0);
        verifyEvaluation("erf(x)+erfc(x)", 2.0, 3.0, 1.0);
        verifyEvaluation("min(3, x)", 2.0, 3.0, 2.0);
        verifyEvaluation("min(y, 5)", 2.0, 3.0, 3.0);
        verifyEvaluation("max(x, y)", 2.0, 3.0, 3.0);
        verifyEvaluation("max(x, -1)", 2.0, 3.0, 2.0);
        verifyEvaluation("abs(x-y)", 2.0, 3.0, 1.0);
        verifyEvaluation("delta(x)+3*delta(y-1.5)", 2.0, 1.5, 3.0);
        verifyEvaluation("step(x-3)+y*step(x)", 2.0, 3.0, 3.0);
        verifyEvaluation("floor(x)", -2.1, 3.0, -3.0);
        verifyEvaluation("ceil(x)", -2.1, 3.0, -2.0);
        verifyEvaluation("select(x, 1.0, y)", 0.3, 2.0, 1.0);
        verifyEvaluation("select(x, 1.0, y)", 0.0, 2.0, 2.0);
        verifyEvaluation("atan2(x, y)", 3.0, 1.5, std::atan(2.0));
        verifyEvaluation("sqrt(x^2)", -2.2, 0.0, 2.2);
        verifyEvaluation("sqrt(x)^2", 2.2, 0.0, 2.2);
        verifyInvalidExpression("1..2");
        verifyInvalidExpression("1*(2+3");
        verifyInvalidExpression("5++4");
        verifyInvalidExpression("1+2)");
        verifyInvalidExpression("cos(2,3)");
        verifyDerivative("x", "1");
        verifyDerivative("x^2+x", "2*x+1");
        verifyDerivative("y^x-x", "log(y)*(y^x)-1");
        verifyDerivative("sin(x)", "cos(x)");
        verifyDerivative("cos(x)", "-sin(x)");
        verifyDerivative("tan(x)", "square(sec(x))");
        verifyDerivative("cot(x)", "-square(csc(x))");
        verifyDerivative("sec(x)", "sec(x)*tan(x)");
        verifyDerivative("csc(x)", "-csc(x)*cot(x)");
        verifyDerivative("exp(2*x)", "2*exp(2*x)");
        verifyDerivative("log(x)", "1/x");
        verifyDerivative("sqrt(x)", "0.5/sqrt(x)");
        verifyDerivative("asin(x)", "1/sqrt(1-x^2)");
        verifyDerivative("acos(x)", "-1/sqrt(1-x^2)");
        verifyDerivative("atan(x)", "1/(1+x^2)");
        verifyDerivative("atan2(2*x,y)", "2*y/(4*x^2+y^2)");
        verifyDerivative("sinh(x)", "cosh(x)");
        verifyDerivative("cosh(x)", "sinh(x)");
        verifyDerivative("tanh(x)", "1/(cosh(x)^2)");
        verifyDerivative("erf(x)", "1.12837916709551*exp(-x^2)");
        verifyDerivative("erfc(x)", "-1.12837916709551*exp(-x^2)");
        verifyDerivative("step(x)*x+step(1-x)*2*x", "step(x)+step(1-x)*2");
        verifyDerivative("recip(x)", "-1/x^2");
        verifyDerivative("square(x)", "2*x");
        verifyDerivative("cube(x)", "3*x^2");
        verifyDerivative("min(x, 2*x)", "step(x-2*x)*2+(1-step(x-2*x))*1");
        verifyDerivative("max(5, x^2)", "(1-step(5-x^2))*2*x");
        verifyDerivative("abs(3*x)", "step(3*x)*3+(1-step(3*x))*-3");
        verifyDerivative("floor(x)+0.5*x*ceil(x)", "0.5*ceil(x)");
        verifyDerivative("select(x, x^2, 3*x)", "select(x, 2*x, 3)");
        testCustomFunction("custom(x, y)/2", "x*y");
        testCustomFunction("custom(x^2, 1)+custom(2, y-1)", "2*x^2+4*(y-1)");
        cout << Parser::parse("x*x").optimize() << endl;
        cout << Parser::parse("x*(x*x)").optimize() << endl;
        cout << Parser::parse("(x*x)*x").optimize() << endl;
        cout << Parser::parse("2*3*x").optimize() << endl;
        cout << Parser::parse("1/(1+x)").optimize() << endl;
        cout << Parser::parse("x^(1/2)").optimize() << endl;
        cout << Parser::parse("log(3*cos(x))^(sqrt(4)-2)").optimize() << endl;
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
