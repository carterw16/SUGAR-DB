"""Voltage control functions for a grid-following IBDG.

  Author(s): Naeem Turner-Bandele
  Created Date: 06-20-2020
  Updated Date: 04-02-2021
  Email: nturnerb@cmu.edu
  Status: Development

  Control the Volt-Var of a grid-following IBDG using a smooth piecewise function, a cubic spline, or an outer loop.

"""

from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scp
import seaborn as sns

sns.set(context="paper", style="whitegrid", font_scale=4)
sns.set_style("ticks")

from lib.scientific.cubic_spline import (create_cubicSpline,
                                         cubicSpline_partials, )
from lib.scientific.polyfix import polyfix


# # # Voltage/Reactive Power Control # # #

class QVSplineControl:

	def __init__(self, V1, V2, V3, V4, Vnom, Qscale):
		self.V1 = V1 * Vnom
		self.V2 = V2 * Vnom
		self.V3 = V3 * Vnom
		self.V4 = V4 * Vnom
		self.V0 = (self.V2 + self.V3) / 2
		self.VL = self.V1 - 0.01 * Vnom
		self.VH = self.V4 + 0.01 * Vnom

		self.Vnom = Vnom
		self.Qscale = Qscale
		self.Qmin = 1 * self.Qscale
		self.Q0 = 0
		self.Qmax = -1 * self.Qscale

		self.Vshift = 0.008 * Vnom
		self.Qshift = 0.20 * self.Qscale
		self.Qmid = 0.10 * self.Qscale

		self.SlopeFactor = 100

		self.Vbounds, self.Qbounds = self.create_bounds()
		self.R1_2, self.R3_4 = self.create_regions()

	def create_bounds(self):
		'''
		Create the bounds for Q-V control function
		:return: Voltage and reactive power bounds
		'''

		# Control the Q Reference setpoint and bounds
		Qmin = self.Qmin
		Q0 = self.Q0
		Qmax = self.Qmax
		Qmid = self.Qmid

		# Create Bounds
		V1_2b = [
			self.V1, (2 * self.V1 + self.Vshift) / 2,
			         (2 * self.V2 - self.Vshift) / 2, self.V2
		]
		Q1_2b = [Qmin, Qmin - Qmid, Q0 + Qmid, Q0]

		# V3_4b = [self.V3, (2 * self.V3 + self.Vshift) / 2, (2 * self.V4 - self.Vshift) / 2, self.V4]
		V3_4b = [
			self.V3, (2 * self.V3 + self.Vshift) / 2,
			         (2 * self.V4 - self.Vshift) / 2, self.V4
		]
		Q3_4b = [Q0, Q0 - Qmid, Qmax + Qmid, Qmax]

		Vbounds = {'V1-2': V1_2b, 'V3-4': V3_4b}
		Qbounds = {'Q1-2': Q1_2b, 'Q3-4': Q3_4b}
		return Vbounds, Qbounds

	def create_regions(self):
		# create cubic spline regions R1, R2, R3, R4
		Q1_2, a1_2, b1_2, c1_2, d1_2 = create_cubicSpline(
			self.Vbounds['V1-2'], self.Qbounds['Q1-2'])
		Q3_4, a3_4, b3_4, c3_4, d3_4 = create_cubicSpline(
			self.Vbounds['V3-4'], self.Qbounds['Q3-4'])

		R1_2 = {'Q': Q1_2, 'a': a1_2, 'b': b1_2, 'c': c1_2, 'd': d1_2}
		R3_4 = {'Q': Q3_4, 'a': a3_4, 'b': b3_4, 'c': c3_4, 'd': d3_4}

		return R1_2, R3_4

	def get_Q(self, V_k=1.0):
		Q0 = self.Q0
		Q1_2 = self.R1_2['Q']
		Q3_4 = self.R3_4['Q']
		Q = 0

		if V_k <= self.V1:
			Q = self.Qmin
		elif self.V1 < V_k <= self.V2:
			Q = 1 * Q1_2(V_k)
		# Q = self.Qscale * Q1_2(V_k)
		elif self.V2 < V_k <= self.V3:
			v2_3 = [self.V2, self.V3]
			Q2_3 = [Q0, Q0]
			m, b, _, _, _ = scp.linregress(v2_3, Q2_3)
			Q = 0
		elif self.V3 < V_k <= self.V4:
			# Q = self.Qscale * Q3_4(V_k)
			Q = 1 * Q3_4(V_k)

		else:
			Q = self.Qmax

		return Q

	def plot_QV(self, V_k=-1.0, save_plot=False):
		rcParams['font.sans-serif'] = "Times New Roman"
		rcParams['font.family'] = "sans-serif"
		rcParams['font.size'] = 36

		font = {'fontname': 'Times New Roman'}

		step = 0.0001
		fig, ax = plt.subplots(figsize=(1280 / 96, 960 / 96), dpi=96)
		ax.grid()

		vl = [self.V1, self.VL]
		Ql = [self.Qmin / 1e6, self.Qmin / 1e6]
		ax.plot(vl, Ql, '-', label='Region 1-2', linewidth=10)

		# # Region 1-2 - V1, Q1  to  V2, Q2 # #
		v1_2 = np.arange(self.V1, self.V2, step)
		Q1_2 = self.R1_2['Q']
		ax.plot(v1_2, Q1_2(v1_2) / 1e6, '-', label='Region 2-3', linewidth=10)

		# # Region 2-3 - Linear Region Between 2 and 3 # #
		v2_3 = [self.V2, self.V3]
		Q2_3 = [self.Q0, self.Q0]
		ax.plot(v2_3, Q2_3, label='Region 3-4', linewidth=10)

		# # Region 3-4 - V3, Q3 to V4, Q4# #
		v3_4 = np.arange(self.V3, self.V4, step)
		Q3_4 = self.R3_4['Q']
		ax.plot(v3_4, Q3_4(v3_4) / 1e6, '-', label='Region 4-5', linewidth=10)

		# Low and High Regions - The upper and lower limits of VL and VH
		vh = [self.V4, self.VH]
		Qh = [self.Qmax / 1e6, self.Qmax / 1e6]
		ax.plot(vh, Qh, label='Region 5-6', linewidth=10)
		# if V_k != -1.0:
		#     if V_k <= self.V1:
		#         ax.plot(V_k, self.Qmin, '*', ms=20.0)
		#     elif self.V1 < V_k <= self.V2:
		#         ax.plot(V_k, Q1_2(V_k), '*', ms=20.0)
		#         # a = self.R2['a']
		#         # b = self.R2['b']
		#         # c = self.R2['c']
		#         # d = self.R2['d']
		#         # Vbt = self.Vbounds['V2'][1]
		#         # F = a[1] * (V_k - Vbt) ** 3 + b[1] * (V_k - Vbt) ** 2 + c[1] * (V_k - Vbt) + d[1]
		#         # ax.plot(V_k, F, 'o', ms=20.0)
		#     elif self.V2 < V_k <= self.V3:
		#         m, b, _, _, _ = scp.linregress(v2_3, Q2_3)
		#         ax.plot(V_k, V_k * m + b, '*', ms=20.0)
		#     elif self.V3 < V_k <= self.V4:
		#         ax.plot(V_k, Q3_4(V_k), '*', ms=20.0)
		#     else:
		#         ax.plot(V_k, self.Qmax, '*', ms=20.0)

		plt.legend(frameon=False)
		plt.ylabel('Reactive Power (MVAr)')
		plt.xlabel('Voltage (V)')
		vticks = [self.VL, self.V1, self.V2, self.V0, self.V3, self.V4, self.VH]
		plt.xticks(
			vticks,
			['$V_L$', '$V_1$', '$V_2$', '$V_{Ref}$', '$V_3$', '$V_4$', '$V_H$'])

		plt.tick_params(width=5, length=15)
		if save_plot:
			plt.savefig('DG_QV_Spline_Control.svg', format='svg', bbox_inches='tight')

		plt.show()

	def control_QV(self, V_k, Vr_k, Vi_k, Q, Qphase='A'):

		flag_Qmax = False
		flag_Qmin = False
		flag_Qlim = False

		FQ = float(0)
		dFQ_dVar = float(0)
		dFQ_dVai = float(0)
		dFQ_dVbr = float(0)
		dFQ_dVbi = float(0)
		dFQ_dVcr = float(0)
		dFQ_dVci = float(0)
		dFQ_dQ = float(1)

		zero = (float(0), float(0))

		if V_k < self.V1:
			flag_Qmin = True
			flag_Qlim = True

		elif self.V1 <= V_k <= self.V2:
			a = self.R1_2['a']
			b = self.R1_2['b']
			c = self.R1_2['c']
			Q1_2 = self.R1_2['Q']
			FQ = Q - Q1_2(V_k)
			if self.Vbounds['V1-2'][0] < V_k <= self.Vbounds['V1-2'][1]:
				p = self.Vbounds['V1-2'][0]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero
			elif self.Vbounds['V1-2'][1] < V_k <= self.Vbounds['V1-2'][2]:
				p = self.Vbounds['V1-2'][1]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero
			else:
				p = self.Vbounds['V1-2'][2]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero

		elif self.V2 < V_k <= self.V3:
			FQ = Q - 0
			dFQ_dQ = 1

		elif self.V3 < V_k <= self.V4:
			a = self.R3_4['a']
			b = self.R3_4['b']
			c = self.R3_4['c']
			d = self.R3_4['d']
			Q3_4 = self.R3_4['Q']
			FQ = Q - Q3_4(V_k)
			if self.Vbounds['V3-4'][0] <= V_k <= self.Vbounds['V3-4'][1]:
				p = self.Vbounds['V3-4'][0]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[0], b[0], c[0], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero
			elif self.Vbounds['V3-4'][1] < V_k <= self.Vbounds['V3-4'][2]:
				p = self.Vbounds['V3-4'][1]
				F = a[1] * (V_k - p) ** 3 + b[1] * (V_k - p) ** 2 + c[1] * (
						V_k - p) + d[1]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[1], b[1], c[1], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero
			else:
				p = self.Vbounds['V3-4'][2]
				dFQ_dVar, dFQ_dVai = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = cubicSpline_partials(
					a[2], b[2], c[2], Vr_k, Vi_k, self.Qscale, p) if Qphase == 'C' else zero

		else:
			flag_Qmax = True
			flag_Qlim = True

		dFQ = {
			'dVr': [dFQ_dVar, dFQ_dVbr, dFQ_dVcr],
			'dVi': [dFQ_dVai, dFQ_dVbi, dFQ_dVci],
			'dQ': dFQ_dQ
		}

		return dFQ, FQ, flag_Qmax, flag_Qmin, flag_Qlim


class QVPatchingControl:
	def __init__(self, V1, V2, V3, V4, Vnom, Qscale):
		self.V1 = V1 * Vnom
		self.V1A = 0.0
		self.V2 = V2 * Vnom
		self.V2A = 0.0
		self.V3 = V3 * Vnom
		self.V3A = 0.0
		self.V4 = V4 * Vnom
		self.V4A = 0.0
		self.V0 = (V2 + V3) / 2
		self.VL = self.V1 - 0.01 * Vnom
		self.VH = self.V4 + 0.01 * Vnom

		self.Vnom = Vnom
		self.Qscale = Qscale
		self.Qmin = 1 * Qscale
		self.Q0 = 0
		self.Q2A = 0.0
		self.Q3A = 0.0
		self.Q4A = 0.0
		self.Qmax = -1 * Qscale

		self.Q_R3 = None
		self.V_R3 = None

		self.Q_R4 = None
		self.V_R4 = None

		self.V_R5 = None

		self.Q_R6 = None
		self.V_R6 = None

		self.Q_R1 = None
		self.V_R1 = None

		self.Q_R7 = None
		self.V_R7 = None

		self.Q_R8 = None
		self.V_R8 = None

		self.Vshift = 0.008 * Vnom
		self.Qshift = 0.20
		self.Qmid = 0.10

		(self.V1A, self.V2A, self.V3A, self.V4A,
		 self.p1, self.p2, self.p3, self.p4,
		 self.Q2A, self.Q3A, self.Q4A) = self.create_bounds()

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
		Q_upper = np.sort(np.linspace(0.995, 1.0, num=51))
		Q_lower = np.sort(np.linspace(0, 0.005, num=51))

		# chosen difference between points
		diff = 0.00050 * self.Vnom

		# Curve Towards Qmax
		y0 = np.array([self.Qmax, self.Qmax])
		p4, xu = self.get_polyval(self.V4, Q_upper, y0, 1, self.Qmax, diff=-diff, sort="ascend")
		yuhat = np.polyval(p4, xu)

		# Curve Towards Qmin
		y0 = np.array([self.Qmin, self.Qmin])
		p1, xl = self.get_polyval(self.V1, Q_upper, y0, 1, self.Qmin, diff=diff, sort="descend")
		ylhat = np.polyval(p1, xl)

		# Upper Curvature Towards Zero
		y0 = np.array([self.Q0, self.Q0])
		p2, xzu = self.get_polyval(self.V2, Q_lower, y0, 1, self.Qmin, diff=-diff, sort="descend")
		yzuhat = np.polyval(p2, xzu)

		# Lower Curvature Away From Zero
		y0 = np.array([self.Q0, self.Q0])
		p3, xzl = self.get_polyval(self.V3, Q_lower, y0, 1, self.Qmax, diff=diff, sort="ascend")
		yzlhat = np.polyval(p3, xzl)

		self.Q_R3 = ylhat
		self.V_R3 = xl

		self.Q_R4 = yzuhat
		self.V_R4 = xzu

		self.V_R5 = np.array([xzu[0], xzl[0]])

		self.Q_R6 = yzlhat
		self.V_R6 = xzl

		self.Q_R1 = np.array([ylhat[0], yzuhat[-1]])
		self.V_R1 = np.array([xl[0], xzu[-1]])

		self.Q_R7 = np.array([yuhat[0], yzlhat[-1]])
		self.V_R7 = np.array([xu[0], xzl[-1]])

		self.Q_R8 = yuhat
		self.V_R8 = xu

		V1A = xl[0]
		V2A = xzu[-1]
		V3A = xzl[-1]
		V4A = xu[0]

		Q2A = yzuhat[-1]
		Q3A = yzlhat[-1]
		Q4A = yuhat[0]

		return V1A, V2A, V3A, V4A, p1, p2, p3, p4, Q2A, Q3A, Q4A

	def get_Q(self, V_k=1.0):

		if V_k <= self.V1:
			Q = self.Qmin
		elif self.V1 < V_k <= self.V1A:
			beta = self.p1
			Q = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V1A < V_k <= self.V2A:
			v1A_2A = [self.V1A, self.V2A]
			Q1A_2A = [1 - self.Q2A, self.Q2A]
			# m, b, _, _, _ = scp.linregress(v1A_2A, Q1A_2A)

			m = (Q1A_2A[1] - Q1A_2A[0]) / (v1A_2A[1] - v1A_2A[0])
			b = (Q1A_2A[0] - (m * self.V2A))

			Q = ((m * V_k) + b)

		elif self.V2A < V_k <= self.V2:
			beta = self.p2
			Q = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V2 < V_k < self.V3:
			Q = 0

		elif self.V3 <= V_k < self.V3A:
			beta = self.p3
			Q = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		elif self.V3A <= V_k < self.V4A:
			v3A_4A = [self.V3A, self.V4A]
			Q3A_4A = [self.Q3A, self.Q4A]

			m = (Q3A_4A[1] - Q3A_4A[0]) / (v3A_4A[1] - v3A_4A[0])
			b = (Q3A_4A[0] - (m * self.V3A))

			Q = ((m * V_k) + b)

		elif self.V4A <= V_k < self.V4:
			beta = self.p4
			Q = ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])

		else:
			Q = self.Qmax

		return Q

	def plot_QV(self, V_k=-1.0, save_plot=False):
		rcParams['font.sans-serif'] = "Times New Roman"
		rcParams['font.family'] = "sans-serif"
		rcParams['font.size'] = 36

		Q1_2 = [np.array([self.Qmin]) / self.Qscale, np.array([self.Qmin]) / self.Qscale]
		V1_2 = [np.array([self.V1]), np.array([self.VL])]
		Q2_3 = self.Q_R3
		V2_3 = self.V_R3
		Q3_4 = self.Q_R1
		V3_4 = self.V_R1

		Q_R4 = self.Q_R4
		V_R4 = self.V_R4
		Q_R5 = np.array([self.Q0, self.Q0])
		V_R5 = self.V_R5
		Q_R6 = self.Q_R6
		V_R6 = self.V_R6
		Q_R7 = self.Q_R7
		V_R7 = self.V_R7
		Q_R8 = self.Q_R8
		V_R8 = self.V_R8
		V_R9 = [np.array([self.V4]), np.array([self.VH])]
		Q_R9 = [np.array([self.Qmax]) / self.Qscale, np.array([self.Qmax]) / self.Qscale]

		fig, ax = plt.subplots(figsize=(1280 / 96, 960 / 96), dpi=96)
		ax.grid()

		region1_2, = ax.plot(V1_2, Q1_2, linewidth=10, label='Region 1-2')
		region2_3, = ax.plot(V2_3, Q2_3 / self.Qscale, linewidth=10, label='Region 2-3')
		region3_4, = ax.plot(V3_4, Q3_4 / self.Qscale, linewidth=10, label='Region 3-4')
		region4_5, = ax.plot(V_R4, Q_R4 / self.Qscale, linewidth=10, label='Region 4-5')
		region5_6, = ax.plot(V_R5, Q_R5 / self.Qscale, linewidth=10, label='Region 5-6')
		region6_7, = ax.plot(V_R6, Q_R6 / self.Qscale, linewidth=10, label='Region 6-7')
		region7_8, = ax.plot(V_R7, Q_R7 / self.Qscale, linewidth=10, label='Region 7-8')
		region8_9, = ax.plot(V_R8, Q_R8 / self.Qscale, linewidth=10, label='Region 8-9')
		region9_10, = ax.plot(V_R9, Q_R9, linewidth=10, label='Region 9-10')

		plt.ylabel('Reactive Power (% Capability)')
		plt.xlabel('Voltage (p.u.)')

		first_legend = plt.legend(handles=[region1_2, region2_3, region3_4, region4_5], loc='upper right',
		                          frameon=False, fontsize='small', markerscale=1.0)
		plt.gca().add_artist(first_legend)

		plt.legend(handles=[region5_6, region6_7, region7_8, region8_9, region9_10], loc='lower left',
		           frameon=False, fontsize='small', markerscale=1.0)

		vticks = [V1_2[1], V1_2[0], (V_R4[0]), V_k,
		          V_R6[0], V_R9[0], V_R9[1]]

		plt.xticks(vticks, ['$V_L$', '$V_1$', '$V_2$', '$V_{Ref}$', '$V_3$', '$V_4$', '$V_H$'])

		plt.tick_params(width=5, length=15)

		if save_plot:
			plt.savefig('DG_QV_Patching_Control.svg', format='svg', bbox_inches='tight')

		plt.show()

	@staticmethod
	def quadratic_partials(Vr, Vi, beta):
		V = np.abs(Vr + Vi * 1j)
		dV_dVr = Vr / V
		dV_dVi = Vi / V

		dF_dV = -(2 * beta[0] * V) - beta[1]

		dF_dVr = dF_dV * dV_dVr
		dF_dVi = dF_dV * dV_dVi

		return dF_dVr, dF_dVi

	@staticmethod
	def linear_partials(Vr, Vi, m):
		V = np.abs(Vr + Vi * 1j)
		dV_dVr = Vr / V
		dV_dVi = Vi / V

		dF_dV = -m

		dF_dVr = dF_dV * dV_dVr
		dF_dVi = dF_dV * dV_dVi

		return dF_dVr, dF_dVi

	def control_QV(self, V_k, Vr_k, Vi_k, Q, Qphase='A'):

		flag_Qmax = False
		flag_Qmin = False
		flag_Qlim = False

		FQ = float(0)
		dFQ_dVar = float(0)
		dFQ_dVai = float(0)
		dFQ_dVbr = float(0)
		dFQ_dVbi = float(0)
		dFQ_dVcr = float(0)
		dFQ_dVci = float(0)
		dFQ_dQ = float(1)

		zero = (float(0), float(0))

		if V_k < self.V1:
			flag_Qmin = True
			flag_Qlim = True

		elif self.V1 <= V_k <= self.V1A:
			beta = self.p1
			FQ = Q - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFQ_dVar, dFQ_dVai = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'C' else zero

		elif self.V1A < V_k <= self.V2A:
			v1A_2A = [self.V1A, self.V2A]
			Q1A_2A = [self.Q_R1[0], self.Q_R1[1]]
			# m, b, _, _, _ = scp.linregress(v1A_2A, Q1A_2A)

			m = (Q1A_2A[1] - Q1A_2A[0]) / (v1A_2A[1] - v1A_2A[0])
			b = (Q1A_2A[0] - (m * self.V1A))

			FQ = Q - ((m * V_k) + b)
			dFQ_dVar, dFQ_dVai = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'C' else zero

		elif self.V2A < V_k <= self.V2:
			beta = self.p2
			FQ = Q - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFQ_dVar, dFQ_dVai = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'C' else zero

		elif self.V2 < V_k < self.V3:
			FQ = Q - 0
			dFQ_dQ = 1

		elif self.V3 <= V_k < self.V3A:
			beta = self.p3
			FQ = Q - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFQ_dVar, dFQ_dVai = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'C' else zero

		elif self.V3A <= V_k < self.V4A:
			v3A_4A = [self.V3A, self.V4A]
			Q3A_4A = [self.Q3A, self.Q4A]

			m = (Q3A_4A[1] - Q3A_4A[0]) / (v3A_4A[1] - v3A_4A[0])
			b = (Q3A_4A[0] - (m * self.V3A))

			FQ = Q - ((m * V_k) + b)
			dFQ_dVar, dFQ_dVai = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'C' else zero

		elif self.V4A <= V_k <= self.V4:
			beta = self.p4
			FQ = Q - ((beta[0] * (V_k ** 2)) + (beta[1] * V_k) + beta[2])
			dFQ_dVar, dFQ_dVai = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'A' else zero
			dFQ_dVbr, dFQ_dVbi = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'B' else zero
			dFQ_dVcr, dFQ_dVci = self.quadratic_partials(Vr_k, Vi_k, beta) if Qphase == 'C' else zero
		else:
			flag_Qmax = True
			flag_Qlim = True

		dFQ = {
			'dVr': [dFQ_dVar, dFQ_dVbr, dFQ_dVcr],
			'dVi': [dFQ_dVai, dFQ_dVbi, dFQ_dVci],
			'dQ': dFQ_dQ
		}

		return dFQ, FQ, flag_Qmax, flag_Qmin, flag_Qlim


class QVLoopControl:

	def __init__(self,
	             Qmax,
	             Qmin,
	             rated_power,
	             Vnom,
	             Va_pu,
	             Vb_pu,
	             Vc_pu,
	             V1,
	             V2,
	             V3,
	             V4):

		self.Qmax = Qmax
		self.Qmin = Qmin
		self.Vbase = Vnom
		self.Va_pu = Va_pu
		self.Vb_pu = Vb_pu
		self.Vc_pu = Vc_pu
		self.rated_power = rated_power

		self.V1 = V1 * self.Vbase
		self.V2 = V2 * self.Vbase
		self.V3 = V3 * self.Vbase
		self.V4 = V4 * self.Vbase

		self.flag_Qmax = False
		self.flag_Qmin = False
		self.flag_Qlim = False

	def linear_partials(self, Vr, Vi, m):
		V = np.abs(Vr + Vi * 1j)
		dV_dVr = Vr / V
		dV_dVi = Vi / V

		dF_dV = -m

		dF_dVr = dF_dV * dV_dVr
		dF_dVi = dF_dV * dV_dVi

		return dF_dVr, dF_dVi

	def get_Q(self, V_k=1.0):
		V_k = V_k / self.Vbase

		if V_k <= self.V1:
			Q = self.Qmin

		elif self.V1 < V_k <= self.V2:
			v1_2 = [self.V1, self.V2]
			Q1_2 = [1, 0]
			m, _, _, _, _ = scp.linregress(v1_2, Q1_2)
			Q = m * V_k * self.Qmin

		elif self.V2 < V_k < self.V3:
			Q = 0

		elif self.V3 <= V_k < self.V4:
			v3_4 = [self.V3, self.V4]
			Q3_4 = [0, -1]
			m, _, _, _, _ = scp.linregress(v3_4, Q3_4)
			Q = m * V_k * self.Qmax

		else:
			Q = self.Qmax

		return Q

	def control_QV(self, V_k, Vr_k, Vi_k, Q, Qphase, flag_Qlim, flag_Qmin, flag_Qmax):

		FQ = float(0)
		dFQ_dVar = float(0)
		dFQ_dVai = float(0)
		dFQ_dVbr = float(0)
		dFQ_dVbi = float(0)
		dFQ_dVcr = float(0)
		dFQ_dVci = float(0)
		dFQ_dQ = float(1)

		zero = (float(0), float(0))

		self.flag_Qmax = flag_Qmax
		self.flag_Qmin = flag_Qmin
		self.flag_Qlim = flag_Qlim

		# V_k = V_k / self.Vbase
		# Vr_k = Vr_k / self.Vbase
		# Vi_k = Vi_k / self.Vbase

		# if V_k < self.V1:
		# 	self.flag_Qmax = True
		# 	self.flag_Qlim = True
		if not self.flag_Qlim:
			if self.V1 <= V_k <= self.V2:
				self.flag_Qmax = False
				self.flag_Qlim = False
				self.flag_Qmin = False

				v1_2 = [self.V1, self.V2]
				Q1_2 = [self.Qmin, 0]
				m = (Q1_2[1] - Q1_2[0]) / (v1_2[1] - v1_2[0])
				b = (Q1_2[0] - (m * self.V1))

				FQ = Q - ((m * V_k) + b)

				dFQ_dVar, dFQ_dVai = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'C' else zero

			elif self.V2 < V_k < self.V3:
				self.flag_Qmax = False
				self.flag_Qlim = False
				self.flag_Qmin = False

				FQ = Q - 0
				dFQ_dQ = 1

			elif self.V3 <= V_k <= self.V4:
				self.flag_Qmax = False
				self.flag_Qlim = False
				self.flag_Qmin = False

				v3_4 = [self.V3, self.V4]
				Q3_4 = [0, self.Qmax]
				m = (Q3_4[1] - Q3_4[0]) / (v3_4[1] - v3_4[0])
				b = (Q3_4[0] - (m * self.V3))
				# m, b, _, _, _ = scp.linregress(v3_4, Q3_4)

				FQ = Q - ((m * V_k) + b)
				dFQ_dVar, dFQ_dVai = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'A' else zero
				dFQ_dVbr, dFQ_dVbi = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'B' else zero
				dFQ_dVcr, dFQ_dVci = self.linear_partials(Vr_k, Vi_k, m) if Qphase == 'C' else zero
		# else:
		# 	self.flag_Qmin = True
		# 	self.flag_Qlim = True

		dFQ = {
			'dVr': [dFQ_dVar, dFQ_dVbr, dFQ_dVcr],
			'dVi': [dFQ_dVai, dFQ_dVbi, dFQ_dVci],
			'dQ': dFQ_dQ
		}

		return dFQ, FQ, self.flag_Qmax, self.flag_Qmin, self.flag_Qlim
