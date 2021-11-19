import constitutive
import numpy as np
import matplotlib.pyplot as plt
import math


##################################################################################################
# Set up properties and model

as4_peek_props = constitutive.Turon06EnrichmentProperties()
# parameters taken from Turon's 06 paper for AS4/PEEK
as4_peek_props.set_parameters(969, 1719, 2.284, 80.0e6, 100.0e6, 1.0e13)
constitutive.EnrichmentPropertiesCollection.singleton().add("as4_peek", as4_peek_props)
props_ref = constitutive.EnrichmentPropertiesCollection.singleton().get_ref(0)
print(props_ref.get().str())

model = constitutive.InfStrainTuron06EnrichmentModel.get()
required_quantities = model.query_required_quantities(props_ref.get())

isvs = constitutive.ISVSetContainer()
model.allocate_isvs(isvs, 1)
state = isvs.get_set(0)

workspace = constitutive.PointwiseEnrichmentModelWorkspace()
model_inputs = constitutive.PointwiseEnrichmentModelInputs()

##################################################################################################
# Perform a simple test to make sure the traction and secant stiffness can be calcualted

solution_jump = np.array([ 0.000007, 0, 0 ])
model_inputs.solution_jump = solution_jump

workspace.initialize(state, model_inputs, props_ref)
model.initialize_model_state(state, model_inputs, props_ref.get(), workspace)
secant_stiffness = model.compute_secant_stiffness_matrix(workspace)
traction         = model.compute_traction(model_inputs.solution_jump, workspace)

print('For solution jump = ', model_inputs.solution_jump);
print('Traction = ', traction)
print(workspace.get_unconverged_state().str())

##################################################################################################
# Plot the contours of the damage variable as a function of the normal and tangential separation

def compute_d(dn, dt):
    """A simple function to compute the damage variable using the model."""
    solution_jump = np.array([ dn, dt / math.sqrt(2.0), dt / math.sqrt(2.0) ])
    model_inputs.solution_jump = solution_jump

    workspace.initialize(state, model_inputs, props_ref)
    model.initialize_model_state(state, model_inputs, props_ref.get(), workspace)
    damage = workspace.get_unconverged_state().get_float_values().get(0)
    
    return damage
    

delta_n, delta_t = np.meshgrid(np.linspace(0, 0.000025, 401),
                               np.linspace(0, 0.000035, 401))

d = np.zeros_like(delta_n)
for i in range(delta_n.shape[0]):                     
    for j in range(delta_n.shape[1]):
        d[i, j] = compute_d(delta_n[i,j], delta_t[i,j])

plt.xlabel('Normal Separation (m)')
plt.ylabel('Tangential Separation (m)')
contour_fill = plt.contourf(delta_n, delta_t, d, levels=30, vmin=0., vmax=1.)
contour_lines = plt.contour(delta_n, delta_t, d, levels=[0.001, 0.1, 0.5, 0.9, 0.999],
                            colors=['white', 'white', 'black', 'black', 'black'], linestyles='dashed')

# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        s = f'{self:.1f}'
        return f'{self:.0f}' if s[-1] == '0' else s

# Recast levels to new class
contour_lines.levels = [nf(val) for val in contour_lines.levels]

# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'd = %r'
else:
    fmt = 'd = %r'

plt.clabel(contour_lines, contour_lines.levels, fmt=fmt, fontsize=10)
plt.colorbar(contour_fill, ticks=[0, 0.5, 1])
plt.show()