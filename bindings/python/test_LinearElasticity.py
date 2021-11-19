import constitutive
import numpy as np
import matplotlib.pyplot as plt

props = constitutive.LinearElasticProperties()
props.set_isotropic_E_nu(1e9, 0.3)
print(props.str())

props2 = constitutive.StateDependentLinearElasticProperties()
props2.set_isotropic_E_nu("1e8 * temperature + 1e8 * x + 1e9", "0.3")
print(props2.str())

c_rel = constitutive.LinearElasticityConstitutiveRelation()
c_rel.set_elastic_properties(props2)

isv = constitutive.LinearElasticityISV()
strain = np.array([0.01, 0, 0, 0, 0, 0])

isv.set_temperature(10)
stiffness_local = np.zeros((6, 6))
props2.update_stiffness_local(isv, [1, 0, 0], stiffness_local)
print(stiffness_local)

sigma_xx = []
T_range = np.linspace(0, 100, 11)
for T in T_range:
    isv.set_temperature(float(T))
    c_rel.update_constitutive_state(isv, [1, 0, 0])
    stress = c_rel.get_stress(isv, strain)
    sigma_xx.append(stress[0][0])

plt.plot(T_range, sigma_xx)
plt.xlabel('Temperature')
plt.ylabel('Stress XX (Pa)')
plt.title('Stress-T response at constant strain')
plt.show()