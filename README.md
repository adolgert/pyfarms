# pyfarms
Simulation of disease spread among farms.

We built this simulation in order to ask how a continuous-time
stochastic simulation of disease spread compares with 
the discrete-time stochastic simulation of disease
spread in the [North American Animal Disease Simulation](http://naadsm.org),
NAADSM.

*Why would you rewrite NAADSM?* There are a host of simulations, not
just of agriculture but of epidemics in people, that write lots of code
based on stories about how a farm works, how people drive around,
how cattle get sick. At the start of the daily time step, the code
walks through that story, flipping a weighted coin whenever there are
options about what might happen next. We think these simulations are
missing a crucial step for translating observation and expert opinion
into a faithful model.

The crucial step is to identify individual behavioral or disease
processes and to mirror them in the code as competing risks for the
next event, the same way it's done in chemical modeling or risk analysis.
This kind of model has many advantages:

 * Can be run in discrete or continuous time.
 * Comparison with other models can be done on a process-by-process basis.
 * Concise description of processes makes these models fast.
 * The rates that come out of a model match those measured in the field.
 * The model specification is distinct from its implementation.
 * Extension of a model means adding processes without disturbing what's there.
 * These models make variance calculations and optimization easier to do.

This code mirrors NAADSM, model for model, which shows that a simulation
as complex as NAADSM can be written in this process-based style.

* Disease model
* Detection and Reporting model
* Quarantine model
* Airborne Spread model
* Indirect and Direct contact models
* Zones are not implemented.

It reads the XML format for NAADSM herd files and scenario files.
(The scenario and herd files from NAADSM are usually missing
a namespace: xmlns:xdf="http://oracc.org/ns/xdf/1.0, which should
be added at the top of the file.)

The Python module for pyfarms depends on the
[PyGSPN](https://github.com/adolgert/PyGSPN) module, also on Dolgert's github.
Installation of this module creates a command `farmspread` which
runs `naadsm.py`. Most of the code is in the file `farms.py`.
