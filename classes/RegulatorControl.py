"""Effective regulator ratio control using patching function.

  Author(s): Naeem Turner-Bandele
  Created Date: 03-31-2021
  Updated Date: 04-02-2021
  Email: nturnerb@cmu.edu
  Status: Development

  Control the effective regulator ratio using a piecewise function where linear and affine regions are linked with
  quadratic patching functions to ensure that the piecewise function is smooth.

"""

from types import SimpleNamespace
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
from lib.scientific.polyfix import polyfix

sns.set(context="paper", style="whitegrid", font_scale=4)
sns.set_style("ticks")


class RegulatorRatioControl:

	def __init__(self, name, Vbase, Vmax, Vmin, max_taps, reg_step):
		self.name = name
		self.Vbase = Vbase
		self.Vmax = Vmax
		self.Vmin = Vmin
		self.max_taps = max_taps / 2
		self.reg_step = reg_step

		self.V1, self.V2, self.V3, self.V4 = self.compute_voltage_bounds()

		# create control regions
		self.aR_R2_R3 = None
		self.V_R2_R3 = None
		self.aR_R4_R5 = None
		self.V_R4_R5 = None
		self.V_R5_R6 = None
		self.aR_R6_R7 = None
		self.V_R6_R7 = None
		self.aR_R3_R4 = None
		self.V_R3_R4 = None
		self.aR_R7_R8 = None
		self.V_R7_R8 = None
		self.aR_R8_R9 = None
		self.V_R8_R9 = None

		(self.V1A, self.V2A, self.V3A, self.V4A,
		 self.p1, self.p2, self.p3, self.p4,
		 self.aR1A, self.aR2A, self.aR3A, self.aR4A) = self.create_bounds()

	def compute_voltage_bounds(self):
		V1 = -self.max_taps * self.reg_step + self.Vmin
		V2 = -self.reg_step + self.Vbase
		V3 = self.reg_step + self.Vbase
		V4 = self.max_taps * self.reg_step + self.Vmax

		return V1, V2, V3, V4

	@staticmethod
	def get_polyval(xlim, yinit, y0, xnom, ynom, diff=-0.005, sort="ascend"):
		xmax = (xlim)
		xmin = (xlim + diff)

		if xmin > xmax:
			x_temp = xmin
			xmin = xmax
			xmax = x_temp

		x = np.sort(np.linspace(xmin, xmax, num=51))

		x = x.reshape((-1, 1))

		if sort == "ascend":
			x = np.sort(xnom * x, axis=0)
		else:
			x = np.sort(xnom * x, axis=0)
			x = x[::-1]

		y = ynom * np.sort(yinit)
		if y0[0] == 0:
			x0 = np.array([x[0], x[2]])
		else:
			x0 = np.array([x[-1], x[-3]])

		p = polyfix(x, y, 2, x0, y0)

		p = p.flatten()

		x = np.reshape(x, (-1, 1))

		return p, x

	def create_bounds(self):
		aR_upper_range = np.sort(np.linspace(0.9985, 1.0, num=51))
		aR_lower_range = np.sort(np.linspace(1.0, 1.005, num=51))
		aR_min_range = np.sort(np.linspace(1.0, 1.01, num=51))
		aR_max_range = aR_upper_range

		aR_max = 1.1
		aR_min = 0.9

		# chosen difference between points
		diff = 0.0050 * self.Vbase
		diff_ascend = 0.0050 * self.Vbase
		# Curve Towards aR min
		y0 = np.array([aR_min, aR_min])
		p1, xl = self.get_polyval(self.V1, aR_min_range, y0, 1, aR_min, diff=diff, sort="descend")
		ylhat = np.polyval(p1, xl)

		# Upper Curvature Towards One
		y0 = np.array([1.0, 1.0])
		p2, xzu = self.get_polyval(self.V2, aR_upper_range, y0, 1, 1, diff=-diff_ascend, sort="ascend")
		yzuhat = np.polyval(p2, xzu)

		# Lower Curvature Away From One
		y0 = np.array([1.0, 1.0])
		p3, xzl = self.get_polyval(self.V3, aR_lower_range, y0, 1, 1, diff=diff, sort="descend")
		yzlhat = np.polyval(p3, xzl)

		# Curve Towards aR max
		y0 = np.array([aR_max, aR_max])
		p4, xu = self.get_polyval(self.V4, aR_max_range, y0, 1, aR_max, diff=-diff_ascend, sort="ascend")
		yuhat = np.polyval(p4, xu)

		self.aR_R2_R3 = ylhat
		self.V_R2_R3 = xl

		self.aR_R4_R5 = yzuhat
		self.V_R4_R5 = xzu

		self.V_R5_R6 = np.array([xzu[-1], xzl[-1]])

		self.aR_R6_R7 = yzlhat
		self.V_R6_R7 = xzl

		self.aR_R3_R4 = np.array([ylhat[0], yzuhat[0]])
		self.V_R3_R4 = np.array([xl[0], xzu[0]])

		self.aR_R7_R8 = np.array([yuhat[0], yzlhat[0]])
		self.V_R7_R8 = np.array([xu[0], xzl[0]])

		self.aR_R8_R9 = yuhat
		self.V_R8_R9 = xu

		# return the additional voltage and aR range points
		V1A = xl[0]
		V2A = xzu[0]
		V3A = xzl[0]
		V4A = xu[0]

		aR1A = ylhat[0]
		aR2A = yzuhat[0]
		aR3A = yzlhat[0]
		aR4A = yuhat[0]

		return V1A, V2A, V3A, V4A, p1, p2, p3, p4, aR1A, aR2A, aR3A, aR4A

	def initialize_aR(self, V_k=120.0):

		# check to see where the load center voltage V_k is within any of the regions below
		# assign an initial aR if the voltage is within that region
		if V_k <= self.V1:
			aR = 0.9
		elif self.V1 < V_k <= self.V1A:
			beta = self.p1
			aR = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V1A < V_k <= self.V2A:
			V_1A_2A = [self.V1A, self.V2A]
			aR_1A_2A = [1 - self.aR2A, self.aR2A]

			m = (aR_1A_2A[1] - aR_1A_2A[0]) / (V_1A_2A[1] - V_1A_2A[0])
			b = (aR_1A_2A[0] - (m * self.V2A))

			aR = ((m * V_k) + b)

		elif self.V2A < V_k <= self.V2:
			beta = self.p2
			aR = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V2 < V_k < self.V3:
			aR = 1.0

		elif self.V3 <= V_k < self.V3A:
			beta = self.p3
			aR = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V3A <= V_k < self.V4A:
			V_3A_4A = [self.V3A, self.V4A]
			aR_3A_4A = [self.Q3A, self.Q4A]

			m = (aR_3A_4A[1] - aR_3A_4A[0]) / (V_3A_4A[1] - V_3A_4A[0])
			b = (aR_3A_4A[0] - (m * self.V3A))

			aR = ((m * V_k) + b)

		elif self.V4A <= V_k < self.V4:
			beta = self.p4
			aR = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		else:
			aR = 1.1

		return aR

	def plot(self, V_k=-1.0, save_plot=False):

		aR_max = 1.1
		aR_min = 0.9

		# Collect the voltage and regulator ratio regions
		aR_R1_R2 = [np.array([aR_min]), np.array([aR_min])]
		V_R1_R2 = [np.array([self.V1]), np.array([self.V1 - 2])]
		aR_R2_R3 = self.aR_R2_R3
		V_R2_R3 = self.V_R2_R3
		aR_R3_R4 = self.aR_R3_R4
		V_R3_R4 = self.V_R3_R4

		aR_R4_R5 = self.aR_R4_R5
		V_R4_R5 = self.V_R4_R5
		aR_R5_R6 = np.array([1.0, 1.0])
		V_R5_R6 = self.V_R5_R6
		aR_R6_R7 = self.aR_R6_R7
		V_R6_R7 = self.V_R6_R7
		aR_R7_R8 = self.aR_R7_R8
		V_R7_R8 = self.V_R7_R8
		aR_R8_R9 = self.aR_R8_R9
		V_R8_R9 = self.V_R8_R9
		V_R9_R10 = [np.array([self.V4]), np.array([self.V4 + 2])]
		aR_R9_R10 = [np.array([aR_max]), np.array([aR_max])]

		# # Plot the regulator control piecewise function with the patching regions # #
		# set the font type and size for the plot
		rcParams['font.sans-serif'] = "Times New Roman"
		rcParams['font.family'] = "sans-serif"
		rcParams['font.size'] = 36

		fig, ax = plt.subplots(figsize=(1280 / 96, 960 / 96), dpi=96)
		ax.grid()

		region1_2, = ax.plot(V_R1_R2, aR_R1_R2, linewidth=10, label='Region 1-2')
		region2_3, = ax.plot(V_R2_R3, aR_R2_R3, linewidth=10, label='Region 2-3')
		region3_4, = ax.plot(V_R3_R4, aR_R3_R4, linewidth=10, label='Region 3-4')
		region4_5, = ax.plot(V_R4_R5, aR_R4_R5, linewidth=10, label='Region 4-5')
		region5_6, = ax.plot(V_R5_R6, aR_R5_R6, linewidth=10, label='Region 5-6')
		region6_7, = ax.plot(V_R6_R7, aR_R6_R7, linewidth=10, label='Region 6-7')
		region7_8, = ax.plot(V_R7_R8, aR_R7_R8, linewidth=10, label='Region 7-8')
		region8_9, = ax.plot(V_R8_R9, aR_R8_R9, linewidth=10, label='Region 8-9')
		region9_10, = ax.plot(V_R9_R10, aR_R9_R10, linewidth=10, label='Region 9-10')

		plt.ylabel("$a_{R}$ (Effective Regulator Ratio)")
		plt.xlabel('Voltage')

		first_legend = plt.legend(handles=[region1_2, region2_3, region3_4, region4_5], loc='lower right',
		                          frameon=False, fontsize='small', markerscale=1.0)
		plt.gca().add_artist(first_legend)

		plt.legend(handles=[region5_6, region6_7, region7_8, region8_9, region9_10], loc='upper left',
		           frameon=False, fontsize='small', markerscale=1.0)

		vticks = [V_R1_R2[0], V_k, V_R9_R10[0]]

		plt.xticks(vticks, ['$V_{L}$', '$V_{center}$', '$V_{H}$'])

		plt.tick_params(width=5, length=15)

		# save the regulator control plot
		if save_plot:
			plt.savefig('Regulator_aR_Control.svg', format='svg', bbox_inches='tight')

		plt.show()

	@staticmethod
	def quadratic_partials(V, beta):
		# V = np.abs(Vr + Vi * 1j)
		# dV_dVr = Vr / V
		# dV_dVi = Vi / V

		dF_dV = -(2 * beta[0] * V) - beta[1]

		# dF_dVr = dF_dV * dV_dVr
		# dF_dVi = dF_dV * dV_dVi

		return dF_dV

	@staticmethod
	def linear_partials(m):
		# V = np.abs(Vr + Vi * 1j)
		# dV_dVr = Vr / V
		# dV_dVi = Vi / V

		dF_dV = -m

		# dF_dVr = dF_dV * dV_dVr
		# dF_dVi = dF_dV * dV_dVi

		return dF_dV

	def compute_linearized_aR(self, V_k, aR, reg_phase):

		# reset / initialize the flags
		flag_aR_max = False
		flag_aR_min = False
		flag_aR_lim = False

		FaR = float(0)
		dFaR_dVa = float(0)
		dFaR_dVb = float(0)
		dFaR_dVc = float(0)
		dFaR_daR = float(1)

		zero = (float(0), float(0))

		if V_k < self.V1:
			flag_aR_min = True
			flag_aR_lim = True

		elif self.V1 <= V_k <= self.V1A:
			beta = self.p1
			FaR = aR - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFaR_dVa = self.quadratic_partials(V_k, beta) if reg_phase == 0 else zero
			dFaR_dVb = self.quadratic_partials(V_k, beta) if reg_phase == 1 else zero
			dFaR_dVc = self.quadratic_partials(V_k, beta) if reg_phase == 2 else zero

		elif self.V1A < V_k <= self.V2A:
			V_1A_2A = [self.V1A, self.V2A]
			aR_1A_2A = [self.aR_R3_R4[0], self.aR_R3_R4[1]]

			m = (aR_1A_2A[1] - aR_1A_2A[0]) / (V_1A_2A[1] - V_1A_2A[0])
			b = (aR_1A_2A[0] - (m * self.V1A))

			FaR = aR - ((m * V_k) + b)
			dFaR_dVa = self.linear_partials(m) if reg_phase == 0 else zero
			dFaR_dVb = self.linear_partials(m) if reg_phase == 1 else zero
			dFaR_dVc = self.linear_partials(m) if reg_phase == 2 else zero

		elif self.V2A < V_k <= self.V2:
			beta = self.p2
			FaR = aR - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFaR_dVa = self.quadratic_partials(V_k, beta) if reg_phase == 0 else zero
			dFaR_dVb = self.quadratic_partials(V_k, beta) if reg_phase == 1 else zero
			dFaR_dVc = self.quadratic_partials(V_k, beta) if reg_phase == 2 else zero

		elif self.V2 < V_k < self.V3:
			FaR = aR - 1.0

		elif self.V3 <= V_k < self.V3A:
			beta = self.p3
			FaR = aR - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFaR_dVa = self.quadratic_partials(V_k, beta) if reg_phase == 0 else zero
			dFaR_dVb = self.quadratic_partials(V_k, beta) if reg_phase == 1 else zero
			dFaR_dVc = self.quadratic_partials(V_k, beta) if reg_phase == 2 else zero

		elif self.V3A <= V_k < self.V4A:
			V_3A_4A = [self.V3A, self.V4A]
			aR_3A_4A = [self.aR3A, self.aR4A]

			m = (aR_3A_4A[1] - aR_3A_4A[0]) / (V_3A_4A[1] - V_3A_4A[0])
			b = (aR_3A_4A[0] - (m * self.V3A))

			FaR = aR - ((m * V_k) + b)
			dFaR_dVa = self.linear_partials(m) if reg_phase == 0 else zero
			dFaR_dVb = self.linear_partials(m) if reg_phase == 1 else zero
			dFaR_dVc = self.linear_partials(m) if reg_phase == 2 else zero

		elif self.V4A <= V_k <= self.V4:
			beta = self.p4
			FaR = aR - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFaR_dVa = self.quadratic_partials(V_k, beta) if reg_phase == 0 else zero
			dFaR_dVb = self.quadratic_partials(V_k, beta) if reg_phase == 1 else zero
			dFaR_dVc = self.quadratic_partials(V_k, beta) if reg_phase == 2 else zero
		else:
			flag_aR_max = True
			flag_aR_lim = True

		dFaR = {
			# 'dVr': [dFaR_dVar, dFaR_dVbr, dFaR_dVcr],
			# 'dVi': [dFaR_dVai, dFaR_dVbi, dFaR_dVci],
			'dV': [dFaR_dVa, dFaR_dVb, dFaR_dVc],
			'daR': dFaR_daR
		}

		dFaR = SimpleNamespace(**dFaR)

		return dFaR, FaR, flag_aR_max, flag_aR_min, flag_aR_lim
