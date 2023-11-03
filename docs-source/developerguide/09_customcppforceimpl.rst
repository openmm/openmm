.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++


.. _customcppforceimpl:

Platform-Independent Code with CustomCPPForceImpl
#################################################

In most cases, Forces are written with platform-specific code that runs as
efficiently as possible on the available hardware.  There are some situations
where that is unnecessarily complicated.  For example, you might write a plugin
that interfaces to an external library through its own public API.  Because all
expensive calculations are done by the external library, there is little
opportunity for any kind of platform-specific optimization.  All you can do is
copy positions over to the external library and copy forces back, hopefully in
a way that is not too inefficient, and preferably with as little code as
possible.  CustomCPPForceImpl is a tool for doing that.

To use it, you write your Force subclass in the usual way, providing whatever
API is appropriate for it.  Then for the corresponding ForceImpl, simply
subclass CustomCPPForceImpl and override :code:`computeForce()` to perform the
calculation using platform-independent C++ code.  Nothing more is required.  It
will automatically work on all platforms.  You do not need to write any GPU
code, provide a Kernel or KernelFactory, define registration functions, or even
create plugin libraries to be loaded dynamically.  The single library containing
the Force and ForceImpl is all you need.

Here is an example of what the code to implement your ForceImpl might look like.
::

    class ExampleForceImpl : public CustomCPPForceImpl {
    public:
        ExampleForceImpl(const ExampleForce& owner) : CustomCPPForceImpl(owner), owner(owner) {
        }

        double computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
            // Compute the forces and energy here.  Store the forces into the
            // vector and return the energy.
        }

        const ExampleForce& getOwner() const {
            return owner;
        }
    private:
        const ExampleForce& owner;
    };
