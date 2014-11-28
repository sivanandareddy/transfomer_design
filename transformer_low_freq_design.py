'''
Created on 19-Nov-2014

@author: sivananda
#Core dimensions here denoted are as per the ferrite EE core. 
#!!!!!Careful while using it for EI core. You need to modify accordingly
#          +--------------------------------------------+
#          |                                            |
#          |    ^      +-------------------+            |
#          |    |      |                   |            |
#          |    |      |                   |            |
#          |    |      |    +--------------+            |
#          |    |      |    |                    ^      |
#          |    |      |    |                    |      |
#          |    |      |    |                    |      |
#          |    |      |    +--------------+ ^   |      |
#          |    +      |                   | |   |      |
#          |    A      |                   | +   +      |
#          |    +      |                   | E   D      |
#          |    |      |                   | +   +      |
#          |    |      |    +--------------+ v   |      |
#          |    |      |    |                    |      |
#          |    |      |    |                    |      |
#          |    |      |    |                    |      |
#          |    |      |    ---------------+     v      |
#          |    |      |                   |            |
#          |    |      |                   |            |
#          |    v      +-------------------+            |
#          |                                            |
#          |                 <-----+F+----->            |
#          |                                            |
#          |           <-------+B+--------->            |
#          +--------------------------------------------+
'''
import math
import cmath

N1 = 312.0;         #Number of primary turns 		<============
N2 = 240.0;           #Number of secondary turns		<============
V1 = 240.0;         #Primary Voltage					<============
V2 = 175.0;           #Secondary Voltage					<============
I2 = 3.0;          #Secondary Current					<============
I1 = I2*V2/V1;      #primary Current
f = 50;             #frequency of operation
C = 2.5*25.4e-3       #Stack height in m					<============
E = 2.0*25.4e-3       #center limb height in m			<============
F = 3.0/2.0*25.4e-3     #half the height of the window <============
B = 5.0/2.0*25.4e-3     #half the height of the core lamination		<============
G = 1.0*25.4e-3		#window width in m				<============
A = 6.0*25.4e-3		#core width in m					<============
pri_gauge_t = 1.219e-3 #primary wire thicknes		<============
pri_gauge_w = 1.219e-3 #primary wire width			<============
sec_gauge_t = 1.422e-3    #secondary winding thickness				<============
sec_gauge_w = 1.422e-3     #secondary winding width					<==========
Ar_phy = C*E        #physical Area of the core
c_s_fac = 0.9			#core stacking factor							<============
Ar = Ar_phy*c_s_fac     #magnetic area of the core Physical Area x Stacking factor
B_max_ach = V1/(4.44*f*Ar*N1)   #Maximum B in the design
win_w = G;          #window width in each lamination
win_h = 2*F;        #window height in each lamination
bob_t_body = 4e-3;  #thickness of the bobin body		<============
bob_t_top  = 3e-3;  #thickness of the bobin top			<============
bob_t_bot  = 3e-3;  #thickness of the bobin bottom		<============
win_w_a = win_w - bob_t_body    # available window width	
win_h_a = win_h - bob_t_bot - bob_t_top # available window height
p_w_axial_sf = 0.8		#primary winding axial stacking factor		<============
s_w_axial_sf = 0.8		#secondary winding axial stacking factor	<============
sheilding_t = 0.5e-3									#sheilding thickness with insulation in m 	<============
u_o = 4*math.pi*1e-7									#permeability of free space
p_w_trans_sf = 0.8									#primary winding transverse stacking factor	<============
s_w_trans_sf = 0.8									#secondary winding transverse stacking factor<============
pri_turns_layer = int((win_h_a - pri_gauge_w)*p_w_axial_sf/pri_gauge_w)	#primary turns per layer
sec_turns_layer = int((win_h_a - sec_gauge_w)*s_w_axial_sf/sec_gauge_w)	#secondary turns per layer
no_pri_layers = N1/pri_turns_layer;				#no. of primary layers calculation
if no_pri_layers > int(no_pri_layers):
    no_pri_layers = int(no_pri_layers)+1;
else:
    no_pri_layers = int(no_pri_layers);
no_sec_layers = N2/sec_turns_layer;				#no. of secondary layers calculation
if no_sec_layers > int(no_sec_layers):
    no_sec_layers = int(no_sec_layers)+1;
else:
    no_sec_layers = int(no_sec_layers);
pri_thick = no_pri_layers * pri_gauge_t/p_w_trans_sf	#thickness of primary
sec_thick = no_sec_layers * sec_gauge_t/s_w_trans_sf	#thickness of secondary
h = pri_thick + sec_thick + sheilding_t;					#total winding thickness
t = h - no_pri_layers*pri_gauge_t - no_sec_layers*sec_gauge_t	#total insulation thickness possible
L_leak = 1.0/3.0*u_o*N1*N1/F/F*(h+2*t)*(F*C+B*(E+2*h))			#L_leakage indutance
MLT_pri = 2.0*(E+2*bob_t_body+pri_thick+C+2*bob_t_body+pri_thick);	#mean turn length of primary
MLT_sec = 2.0*(E+2*bob_t_body+2*pri_thick+2*sheilding_t+sec_thick+ C+2*bob_t_body+2*pri_thick+2*sheilding_t+sec_thick);	#mean turn length of secondary
app_pri_len = MLT_pri*N1;	#approximate winding length of primary
app_sec_len = MLT_sec*N2;	#approximate winding length of secondary
con_ar_pri = math.pi*pri_gauge_t*pri_gauge_t/4	#conductor crossection area of primary
#con_ar_sec = math.pi*sec_gauge_t*sec_gauge_t/2;				#conductor cross-section area of secondary
#con_ar_sec = sec_gauge_t*sec_gauge_w;			#conductor cross-section area of secondary
con_ar_sec = math.pi*sec_gauge_t*sec_gauge_t/4;				#conductor cross-section area of secondary
cur_dense_pri = I1/con_ar_pri*1e-6					#primary current density 
cur_dense_sec = I2/con_ar_sec*1e-6					#secondary current density
ro = 1.7e-8;												#Resistivity of copper 1.7x10^8
R_pri = ro*app_pri_len/con_ar_pri					#Resistance of copper for primary
R_sec = ro*app_sec_len/con_ar_sec					#Resistance of copper for secondary
R_sec_1 = R_sec * N1/N2 * N1/N2						#Resistance of copper for secondary when refered to primary
P_cu_pri = I1*I1*(R_pri);								#primary copper loses
P_cu_sec = I2*I2*(R_sec);								#Secondary copper losses
P_cu_tot = P_cu_pri + P_cu_sec;						#Total copper losses
Z_source = R_sec_1 + R_pri + 2*math.pi*f*1j*L_leak;	#Thevenins impedence of source
R_load_pri = V2/I2*N1/N2*N1/N2						#output load when refered to primary
V_out = R_load_pri/(Z_source+R_load_pri)*V1;		#output voltage when dropped
V_out_abs = abs(V_out);									#Abs(output voltage)
V_out_phase_deg = cmath.phase(V_out)*180.0/math.pi; #angle(output voltage)
reg = (V1 - V_out_abs)*100/V1							#Regulation
vol_core = (A*2*B - 2*2*F*G)*C;						#Volume of the core
density = 7650;											#Density of CRGO core in Kg/cu.m
wt_core = vol_core * density;							#Weight of the core in Kg
P_fe = wt_core * 1.11; 									#Iron loss (P_fe per Kg = 1.11W)
surf_A_fe = 2*(2*B+A)*C + (B*A) - (E+2*G)*2*F	#Surface area of Iron 
surf_A_cu = 2*F*(E+2*h+2*bob_t_body)				#Surface area of copper
heat_dis_fe = P_fe/surf_A_fe;							#heat dissipation in Iron
heat_dis_cu = P_cu_tot/surf_A_cu; 					#heat dissipation in Cu
density_cu = 8960;
vol_pri_cu = con_ar_pri*app_pri_len
wt_pri_cu = vol_pri_cu*density_cu;
vol_sec_cu =  con_ar_sec*app_sec_len;
wt_sec_cu = vol_sec_cu*density_cu;

print ("\n #============================================================#")
print ("Remarks: Stack height 2.5 inches, 43No , Transformer - 4;as designed by supplier")
print ("Core stacking factor = "+str(c_s_fac))
print ("primary voltage V1 = "+str(V1)+" V");
print ("Secondary voltage V2 = " + str(V2)+" V")
print ("Secondary Voltage calculated from turns ratio =" +str(V1/N1*N2)+" V")
print ("Primary Current I1 = " +str(I1)+" A")
print ("Secondary Current I2 = "+str(I2)+" A")
print ("primary turns per layer = "+str(pri_turns_layer)+" turns per layer")
print ("Primary winding axial stacking factor = " +str(p_w_axial_sf))
print ("primary turns = "+str(N1)+" turns")
print ("primary layer = "+str(no_pri_layers)+ " layers")
print ("Secondary turns per layer = "+str(sec_turns_layer)+ " turns per layer")
print ("Secondary winding axial stacking factor = " +str(s_w_axial_sf))
print ("Seconcary turns = "+str(N2)+" turns")
print ("Secondary layer = "+str(no_sec_layers)+ " layers")
print ("Primary winding transverse stacking factor = " +str(p_w_trans_sf))
print ("Secondary winding transverse stacking factor = " +str(s_w_trans_sf))
print ("Primary winding thickness = "+str(pri_thick*1000)+" mm")
print ("Secondary winding thickness = "+str(sec_thick*1000)+" mm")
print ("B_max = "+str(B_max_ach)+" T");
print ("MLT primary = "+str(MLT_pri)+" m")
print ("MLT secondary = "+str(MLT_sec)+" m")
print ("Resistance primary = "+str(R_pri)+" Ohms")
print ("Resistance secondary = "+str(R_sec)+" Ohms")
print ("Resistance secondary_ref = "+str(R_sec_1)+" Ohms")
print ("Leakage inductance = "+str(L_leak)+" H")
print ("Copper_losses = "+str(P_cu_tot)+" W")
print ("Heat Dissipation in Cu = "+str(heat_dis_cu)+" W/sq.m");
print ("primary copper losses = "+str(P_cu_pri)+" W")
print ("secondary copper losses = "+str(P_cu_sec)+" W")
print ("Core losses = "+str(P_fe)+" W")
print ("Heat Dissipation in iron = "+str(heat_dis_fe)+" W/sq.m");
print ("Regulation = "+str(reg)+" %")
print ("Primary current density = "+str(cur_dense_pri)+" A/sq.mm")
print ("Secondary current density = "+str(cur_dense_sec)+" A/sq.mm")
print ("Primary winding Length = "+str(app_pri_len)+" m")
print ("Secondary winding Length = "+str(app_sec_len)+" m")
print ("total winding thickness = "+str(h)+" m")
print ("Window width = "+str(G)+" m")
print ("primary gauge = "+str(18)+" SWG")
print ("secondary gauge = "+str(17)+" SWG")
#print ("secondary gauge = 3 strips of size 7mm x 2.4mm")
print ("weight of core = "+str(wt_core)+" kg");
print ("weight of primary copper = "+str(wt_pri_cu)+" kg")
print ("weight of secondary copper = "+str(wt_sec_cu)+" kg");
print ("weight of total copper = "+str(wt_pri_cu+wt_sec_cu)+" kg");
print ("weight of transformer = "+str(wt_pri_cu+wt_sec_cu+wt_core)+" kg")
print ("#============================================================#")
