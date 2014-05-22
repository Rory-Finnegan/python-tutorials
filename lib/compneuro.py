"""
Most of this module is based of the hodgkin_huxley example found here.
http://briansimulator.org/demo/
"""

from numpy.random import choice
from brian import Equations, NeuronGroup, EmpiricalThreshold, Connection, StateMonitor, run, plot, show, randn
from brian import umetre, ufarad, cm, mV, siemens, msiemens, nS, ms, msecond

def hodgkin_huxley(duration=100, num=10, percent_excited=0.7, sample=1):
    """
    The hodgkin_huxley function takes the following parameters:

      duration - is the timeperiod to model for measured in milliseconds.
      num - an integer value represeting the number of neurons you want to model.
      percent_excited - a float between 0 and 1 representing the percentage of excited neurons.
      sample - gives the number of random neurons you would like plotted (default is 1)
    """
    # Assert that we are getting valid input
    assert(percent_excited >= 0 and percent_excited <= 1.0)
    assert(duration > 0)
    assert(num > 0)
    assert(sample > 0)
    assert(num >= sample)

    # Constants used in the modeling equation
    area = 20000*umetre**2
    Cm = (1*ufarad*cm**-2)*area
    gl = (5e-5*siemens*cm**-2)*area
    El = -60*mV
    EK = -90*mV
    ENa = 50*mV
    g_na = (100*msiemens*cm**-2)*area
    g_kd = (30*msiemens*cm**-2)*area
    VT = -63*mV

    # Time constants
    taue = 5*ms       # excitatory
    taui = 10*ms      # inhibitory
    
    # Reversal potentials
    Ee = 0*mV         # excitatory
    Ei = -80*mV       # inhibitory
    
    # Synaptic weights
    we = 6*nS         # excitatory
    wi = 67*nS        # inhibitory

    # The model equations
    eqs = Equations(
        '''
        dv/dt = (gl*(El-v)+ge*(Ee-v)+gi*(Ei-v)-g_na*(m*m*m)*h*(v-ENa)-g_kd*(n*n*n*n)*(v-EK))/Cm : volt
        dm/dt = alpham*(1-m)-betam*m : 1
        dn/dt = alphan*(1-n)-betan*n : 1
        dh/dt = alphah*(1-h)-betah*h : 1
        dge/dt = -ge*(1./taue) : siemens
        dgi/dt = -gi*(1./taui) : siemens
        alpham = 0.32*(mV**-1)*(13*mV-v+VT)/(exp((13*mV-v+VT)/(4*mV))-1.)/ms : Hz
        betam = 0.28*(mV**-1)*(v-VT-40*mV)/(exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
        alphah = 0.128*exp((17*mV-v+VT)/(18*mV))/ms : Hz
        betah = 4./(1+exp((40*mV-v+VT)/(5*mV)))/ms : Hz
        alphan = 0.032*(mV**-1)*(15*mV-v+VT)/(exp((15*mV-v+VT)/(5*mV))-1.)/ms : Hz
        betan = .5*exp((10*mV-v+VT)/(40*mV))/ms : Hz
        '''
    )

    # Build the neuron group
    neurons = NeuronGroup(
        num,
        model=eqs,
        threshold=EmpiricalThreshold(threshold=-20*mV,refractory=3*ms),
        implicit=True,
        freeze=True
    )

    num_excited = int(num * percent_excited)
    num_inhibited = num - num_excited
    excited = neurons.subgroup(num_excited)
    inhibited = neurons.subgroup(num_inhibited)
    excited_conn = Connection(excited, neurons, 'ge', weight=we, sparseness=0.02)
    inhibited_conn = Connection(inhibited, neurons, 'gi', weight=wi, sparseness=0.02)
    
    # Initialization
    neurons.v = El+(randn(num)*5-5)*mV
    neurons.ge = (randn(num)*1.5+4)*10.*nS
    neurons.gi = (randn(num)*12+20)*10.*nS
    
    # Record a few trace and run
    recorded = choice(num, sample)
    trace = StateMonitor(neurons, 'v', record=recorded)
    run(duration * msecond)

    for i in recorded:
        plot(trace.times/ms, trace[i]/mV)
        
    show()

if __name__ == '__main__':
    hodgkin_huxley(percent_excited=0.6, sample=5)

