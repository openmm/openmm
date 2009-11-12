/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

/**
 * This file contains the routine for evaluating a custom expression.
 */

static __constant__ float globalParams[8];

texture<float4, 1, cudaReadModeElementType> texRef0;
texture<float4, 1, cudaReadModeElementType> texRef1;
texture<float4, 1, cudaReadModeElementType> texRef2;
texture<float4, 1, cudaReadModeElementType> texRef3;

#define STACK(y) stack[(y)*blockDim.x+threadIdx.x]
#define VARIABLE(y) variables[(y)*blockDim.x+threadIdx.x]

template<int SIZE>
__device__ float kEvaluateExpression_kernel(Expression<SIZE>* expression, float* stack, float* variables)
{
    int stackPointer = -1;
    for (int i = 0; i < expression->length; i++)
    {
        int op = expression->op[i];
        if (op < MULTIPLY) {
            STACK(++stackPointer) = VARIABLE(op-VARIABLE0);
        }
        else if (op < NEGATE) {
            if (op < MULTIPLY_CONSTANT) {
                if (op == MULTIPLY) {
                    float temp = STACK(stackPointer);
                    STACK(--stackPointer) *= temp;
                }
                else if (op == DIVIDE) {
                    float temp = STACK(stackPointer);
                    STACK(stackPointer) = temp/STACK(--stackPointer);
                }
                else if (op == ADD) {
                    float temp = STACK(stackPointer);
                    STACK(--stackPointer) += temp;
                }
                else if (op == SUBTRACT) {
                    float temp = STACK(stackPointer);
                    STACK(stackPointer) = temp-STACK(--stackPointer);
                }
                else /*if (op == POWER)*/ {
                    float temp = STACK(stackPointer);
                    STACK(stackPointer) = pow(temp, STACK(--stackPointer));
                }
            }
            else if (op < GLOBAL) {
                if (op == MULTIPLY_CONSTANT) {
                    STACK(stackPointer) *= expression->arg[i];
                }
                else if (op == POWER_CONSTANT) {
                    STACK(stackPointer) = pow(STACK(stackPointer), expression->arg[i]);
                }
                else /*if (op == ADD_CONSTANT)*/ {
                    STACK(stackPointer) += expression->arg[i];
                }
            }
            else {
                if (op == GLOBAL) {
                    STACK(++stackPointer) = globalParams[(int) expression->arg[i]];
                }
                else if (op == CONSTANT) {
                    STACK(++stackPointer) = expression->arg[i];
                }
                else /*if (op == CUSTOM || op == CUSTOM_DERIV)*/ {
                    int function = (int) expression->arg[i];
                    float x = STACK(stackPointer);
                    float4 params = cSim.pTabulatedFunctionParams[function];
                    if (x < params.x || x > params.y)
                        STACK(stackPointer) = 0.0f;
                    else
                    {
                        int index = floor((x-params.x)*params.z);
                        float4 coeff;
                        if (function == 0)
                            coeff = tex1Dfetch(texRef0, index);
                        else if (function == 1)
                            coeff = tex1Dfetch(texRef1, index);
                        else if (function == 2)
                            coeff = tex1Dfetch(texRef2, index);
                        else
                            coeff = tex1Dfetch(texRef3, index);
                        x = (x-params.x)*params.z-index;
                        if (op == CUSTOM)
                            STACK(stackPointer) = coeff.x+x*(coeff.y+x*(coeff.z+x*coeff.w));
                        else
                            STACK(stackPointer) = (coeff.y+x*(2.0f*coeff.z+x*3.0f*coeff.w))*params.z;
                    }
                }
            }
        }
        else {
            if (op < SIN) {
                if (op == NEGATE) {
                    STACK(stackPointer) *= -1.0f;
                }
                else if (op == RECIPROCAL) {
                    STACK(stackPointer) = 1.0f/STACK(stackPointer);
                }
                else if (op == SQRT) {
                    STACK(stackPointer) = sqrt(STACK(stackPointer));
                }
                else if (op == EXP) {
                    STACK(stackPointer) = exp(STACK(stackPointer));
                }
                else if (op == LOG) {
                    STACK(stackPointer) = log(STACK(stackPointer));
                }
                else if (op == SQUARE) {
                    float temp = STACK(stackPointer);
                    STACK(stackPointer) *= temp;
                }
                else /*if (op == CUBE)*/ {
                    float temp = STACK(stackPointer);
                    STACK(stackPointer) *= temp*temp;
                }
            }
            else {
                if (op == SIN) {
                    STACK(stackPointer) = sin(STACK(stackPointer));
                }
                else if (op == COS) {
                    STACK(stackPointer) = cos(STACK(stackPointer));
                }
                else if (op == SEC) {
                    STACK(stackPointer) = 1.0f/cos(STACK(stackPointer));
                }
                else if (op == CSC) {
                    STACK(stackPointer) = 1.0f/sin(STACK(stackPointer));
                }
                else if (op == TAN) {
                    STACK(stackPointer) = tan(STACK(stackPointer));
                }
                else if (op == COT) {
                    STACK(stackPointer) = 1.0f/tan(STACK(stackPointer));
                }
                else if (op == ASIN) {
                    STACK(stackPointer) = asin(STACK(stackPointer));
                }
                else if (op == ACOS) {
                    STACK(stackPointer) = acos(STACK(stackPointer));
                }
                else if (op == ATAN) {
                    STACK(stackPointer) = atan(STACK(stackPointer));
                }
                else if (op == SINH) {
                    STACK(stackPointer) = sinh(STACK(stackPointer));
                }
                else if (op == COSH) {
                    STACK(stackPointer) = cosh(STACK(stackPointer));
                }
                else /*if (op == TANH)*/ {
                    STACK(stackPointer) = tanh(STACK(stackPointer));
                }
            }
        }
    }
    return STACK(stackPointer);
}
