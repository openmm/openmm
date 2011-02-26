%pythonprepend %{
try: args=stripUnits(args)
except UnboundLocalError: pass
%}

