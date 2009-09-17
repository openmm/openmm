/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "OpenCLPlatform.h"
#include "OpenCLStreamFactory.h"
#include "OpenCLStreamImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"

using namespace OpenMM;

StreamImpl* OpenCLStreamFactory::createStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, ContextImpl& context) const {
    OpenCLPlatform::PlatformData& data = *static_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData());
    if (name == "particlePositions") {
        float padding[] = {100000.0f, 100000.0f, 100000.0f, 0.2f};
        return new OpenCLStreamImpl<cl_float4>(name, size, type, platform, &data.context->getPosq(), 4, *data.context);
    }
    if (name == "particleVelocities") {
        float padding[] = {0.0f, 0.0f, 0.0f, 0.0f};
        return new OpenCLStreamImpl<cl_float4>(name, size, type, platform, &data.context->getVelm(), 4, *data.context);
    }
    if (name == "particleForces") {
        float padding[] = {0.0f, 0.0f, 0.0f, 0.0f};
        return new OpenCLStreamImpl<cl_float4>(name, size, type, platform, &data.context->getForce(), 4, *data.context);
    }
    switch (type) {
    case Stream::Float:
    case Stream::Double:
        return new OpenCLStreamImpl<cl_float>(name, size, type, platform, *data.context);
    case Stream::Float2:
    case Stream::Double2:
        return new OpenCLStreamImpl<cl_float>(name, size, type, platform, *data.context);
    case Stream::Float3:
    case Stream::Double3:
        return new OpenCLStreamImpl<cl_float>(name, size, type, platform, *data.context);
    case Stream::Float4:
    case Stream::Double4:
        return new OpenCLStreamImpl<cl_float>(name, size, type, platform, *data.context);
    case Stream::Integer:
        return new OpenCLStreamImpl<cl_int>(name, size, type, platform, *data.context);
    case Stream::Integer2:
        return new OpenCLStreamImpl<cl_int>(name, size, type, platform, *data.context);
    case Stream::Integer3:
        return new OpenCLStreamImpl<cl_int>(name, size, type, platform, *data.context);
    case Stream::Integer4:
        return new OpenCLStreamImpl<cl_int>(name, size, type, platform, *data.context);
    }
    throw OpenMMException("Tried to create a Stream with an illegal DataType.");
}
