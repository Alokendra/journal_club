# Active Volume Regulation of Adhered Cells
### Ram Adar, Sam Saffran
### PNAS 2020


## Introduction

  - Motivated by recent experiments showing volume of adhered cells is *reduced* as their basal area is *increased*.
  - Factors considered in the model are:
	- Mechanical equilibrium with extracellular environment.
	- Electrostatic equilibrium through active ion transport.
  - Model relates cell volume with basal area and osmotic pressure.
  
Note:
1) During spreading, the cell volume decreases by several thousand cubic micrometers, corresponding to large pressure changes of the order of megapascals

2) These factors are manifested through two mesoscale, measurable quantities.


## Cell Volume Control

<img src="presentation/images/Cell_Volume_Control.jpg" alt="Cell Volume Control" style="width: 70%, height: 70%">

<div class="references">
Ginzberg, Miriam B., Ran Kafri, and Marc Kirschner. "On being the right (cell) size." Science 348.6236 (2015): 1245075.
</div>

<p style="font-size:70%;"> Cell size is not the result of physical constraint but is adaptively regulated.</p>

Note:
Sizes of different human cell types. Cells are shown to scale. Pancreatic beta cells (insulin and DNA stained) (76), hepatocytes (β-catenin and DNA stained) (77), keratinocyes from oral tissue (78), fibroblasts (79), adipocytes from subcutaneous tissue (80)


## Cell Volume and Cell Cycle

<img src="presentation/images/Cell_Volume_Cell_Cycle.jpg" width=60% height=60%>
<div class="references">
Ginzberg, Miriam B., Ran Kafri, and Marc Kirschner. Science 348.6236 (2015): 1245075.
</div>

<div class="figurecaption">
A critical cell mass is needed for S-phase entry.
</div>

Note:
In populations of proliferating cells, size uniformity may be ensured by linking the processes of growth and cell cycle progression. One way this can be accomplished is by restricting progress through a particular cell cycle stage (for example, the G1/S transition) to cells that have reached a specific “target” size.


## Functional Significance

<img src="presentation/images/Cell_Volume_Functional_Significance.jpeg" width=60% height=60%>
<div class="references">
Lang, Florian. Journal of the American college of nutrition 26.sup5
(2007): 613S-623S.
</div>

<div class="figurecaption">
Cell volume and volume sensitive functions participate in a wide variety of physiological mechanisms.
</div>

Note:
Even at constant extracellular osmolarity, the cell volume can change due to change in extracellular fluid composition. Entry of K+, bicarbonate or organic anions can cause in cell swelling. It can also change due to tranport of osmotically active substances.


# Theories of Cell Volume Regulation



## Liquid Drop Model

<img src="presentation/images/Cell_Volume_Laplace_Law.png" width=60% height=60%>
<div class="references">
Clark, Andrew G., and Ewa Paluch. Springer, Berlin, Heidelberg, 2011. 31-73.
(2007): 613S-623S.
</div>

In most cells surface tension $ \gamma << T $ and cortical tension $ T $ generates hydrostatic pressure.

Note:
Surface tension, intracellular pressure and implications for the regulation of the cell
volume. (a) In most cells, a cortical network of cross-linked actin filaments lies under the plasma
membrane. The surface tension of the cell is a combination of the tension of the plasma membrane,
g, and the cortical tension, T. If the cortex is tightly attached to the membrane, the total surface
tension is simply $\gamma + T$. In most cells, $\gamma$ is orders of magnitude smaller than T, and g can thus be
neglected. (b) Cortical tension generates a hydrostatic pressure in the cytoplasm. The law of
Laplace relates the difference between intracellular pressure, P in , and external pressure, P out , to the
cell radius, R, and the cell surface tension, T


### Osmotic Pressure, Tension, Volume

* Hydrostatic Pressure Balance:
<div class="tooltip">$P_{in} - P_{out} = \frac{2T}{R_c}$
   <span class="tooltiptext">Rc = critical radius, T = cortical tension</span>
</div>

* Osmotic pressure Balance
<div class="tooltip">$\Pi_{in} - \Pi_{out} = \frac{n_{osm}R\Theta}{V_c} - \Pi_{out}$
   <span class="tooltiptext">nₒₛₘ = osmotic pressure, R = gas constant, V = cell volume</span>
</div>


* At equilibrium $P_{in} - P_{out} = \Pi_{in} - \Pi_{out}$ which gives
$$V_c = \frac{n_{osm}R\Theta}{\Pi_{out} + \frac{2T}{R_c}}$$

* Cortical tension $\Pi_{out} \approx 10^5 Pa >> 2T/R_c \approx 10^3 Pa$ which gives

$$V_c \approx \frac{n_{osm}R\Theta}{\Pi_{out}}$$



## Ion Transport

<img src="presentation/images/Ion_Transport.jpg" width=70% height=70%>
<div class="references">
Jiang, Hongyuan, and Sean X. Sun. Biophysical journal 105.3 (2013): 609-619.
</div>

Note:
(a) A minimal model of cell volume and pressure regulation.
We consider a spherical cell enclosed by an actomyosin cortex and the cell
membrane. Embedded in the membrane are several families of passive MS
ion channels and active ion pumps. The MS channels and active ion pumps
can change the internal ion concentration, c in , and the intneral osmotic pres-
sure, P in , leading to changes in water flux across the membrane. Net flow of
water leads to cell volume changes. 
(b) Opening probability versus cortical
stress for an MS channel. The red curve shows the simplified piecewise
linear function used in Eq. 2.
(c) The flux of ions transported by ion trans-
porters as a function of osmotic pressure difference, DP. This flux is
modeled by Eq. 3. (d) The steady-state phase diagram of the cell as a func-
tion of model parameters DP c and s c . The model predicts two regimes. In
the static regime, there are no ion fluxes at steady-state cell size. All of the
channels and pumps are inactive. In the pump-and-leak regime, the influx
and efflux of ions balance, and the cell maintains a steady-state size.


### Cell Volume and Water Flux

* Flux of water:
<div class="tooltip">$J_{water} = -\alpha\Delta\Psi = -\alpha(\Delta P -\Delta\Pi)$
   <span class="tooltiptext">$\alpha$ is a rate constant, $\Delta\Psi$ is chemical potential difference</span>
</div>

* Relation with cell volume and surface area:
<div class="tooltip">$J_{water} = -\alpha(\Delta P -\Delta\Pi) = \frac{1}{A}\frac{dV}{dt}$
   <span class="tooltiptext">For spherical cells $\frac{1}{A}\frac{dV}{dt} = \frac{dr}{dt}$</span>
</div>


### Flux of ions

* Flux across N channels:
<div class="tooltip">$$J_1 \propto NP_{open}\frac{\Delta c}{h_0} = \frac{NP_{open}\Delta\Pi}{RTh_\theta}$$
   <span class="tooltiptext">$\Delta c/h_0$ = Concentration gradient, $P_{open}$ = Membrane opening probability</span>
</div>

* Flux across pumps:
<div class="tooltip">
$$J_2 = -\gamma' RT\ln\frac{c_{in}}{c_{out}} - \Delta G_a $$
   <span class="tooltiptext">$\gamma'$ = Another rate, $\Delta G_a$ is the activation energy needed to pump the ions</span>
</div>


* Considering small concentration gradient $(c_{in} - c_{out})/c_{in} << 1$
<div class="tooltip">
$$J_2 \approx RT\frac{\Pi_{in} - \Pi_{out}}{\Pi_{out}} - \Delta G_a$$
   <span class="tooltiptext">$\Pi = cRT$</span>
</div>

* Net rate of ion transport

<div class="tooltip">
$$\frac{dn}{dt} = A(J_1 + J_2) = A((\beta + \gamma)\Delta\Pi - \gamma\Delta\Pi_a)$$
   <span class="tooltiptext">$\Delta\Pi_a = \frac{\Pi_{out}\Delta G_a}{RT}$</span>
</div>



# Results from Cited Experimental Papers


## Cell Volume and Stiffness

<img src="presentation/images/Guo_Volume_vs_Area.jpg" width=50% height=50%>
<div class="references">
Guo, Ming, et al. Proceedings of the National Academy of Sciences 114.41 (2017): E8618-E8627.
</div>

Note:
This paper considers two time scales of cell volume change - one through change of number of proteins
and cell growth and the other a more rapid change caused by transport of cells in restricted spaces.
(A) Shown are 3D images of A7 cells on micropatterned islands of different sizes on glass. Cells are labeled with cell tracker green.
(Scale bars, 20 μm.) (B) Cell volume decreases with increasing cell spread area on glass. (C) Cell volume plotted as a function of the projected area, for cells on
substrates with different stiffnesses (gray circles; n > 200), cells on a glass substrate but with different available spread area (blue squares; n > 200), and a
dynamically spreading cell (red crosses; n = 3). (D) Variation of cell spread area and volume as a single cell dynamically attaches on a stiff substrate (n = 3).
(E) Schematic illustration of cell volume decrease through water efflux, as cells spread out or are osmotically compressed.


### Bulk Modulus

<img src="presentation/images/Guo_Fig1.png" width=45% height=45%>
<div class="tooltip">$B = -VdP/dV = Nk_BT/(V - V_{min})^2$
   <span class="tooltiptext">$N$ = Number of intracellular osmolytes, $V_{min}$ = Min volume under extreme compression</span>
</div>


### Nuclear Volume

<img src="presentation/images/Guo_Fig2.png" width=70% height=70%>
<div class="references">
Guo, Ming, et al. Proceedings of the National Academy of Sciences 114.41 (2017): E8618-E8627.
</div>


## Controlling Cell Volume via Substrate

<img src="presentation/images/Xie_Cell_Volume_Stiffness.jpg" width=90% height=90%>
<div class="references">
Xie, Kenan, Yuehua Yang, and Hongyuan Jiang.  Biophysical journal 114.3 (2018): 675-687.
</div>

Note:
From the above experiments, we can see that the bigger the
spread area is, the smaller the cell volume and cell height
become. Strikingly, when all the data of living 3T3 cells
are plotted in the same figure, we find that both the cell vol-
ume and cell height decay exponentially with the spread
area (Fig. 3, a and b) no matter what method is used to
decrease cell volume. Therefore, when cells adhere to
substrate, their volume can be determined by the spread
area. It should be noted that the cell volume of suspension
cells is also on the exponential curve.


<img src="presentation/images/Xie_Cell_Volume_Stiffness2.jpg" width=70% height=70%>
<div class="references">
Xie, Kenan, Yuehua Yang, and Hongyuan Jiang.  Biophysical journal 114.3 (2018): 675-687.
</div>



# Article Discussion


## Model Description

<img src="presentation/images/Model_Description.jpg" width=70% height=70%>
<div class="references">
Adar, Ram M., and Samuel A. Safran. Proceedings of the National Academy of Sciences 117.11 (2020): 5604-5609.
</div>

* Model considers proteins and their *excluded* volume.
* Nuclear volume considered to be *fixed*.
* Extracellular medium considered to be *buffered* with a single species.


* Consider cell has $N$ proteins and nucleus occupies a fraction $\alpha$ of the cell volume $V$.

* Effective concentration of proteins $n_p$ considering protein excluded volume and nuclear volume is

<div class="tooltip">$n_p = \frac{N}{(1 - \alpha)V - Nv} = \frac{N_p}{V - V_m}$
   <span class="tooltiptext">$N_p = N/(1 - \alpha)$ Effective protein volume , $V_m = N_pv$ minimal volume</span>
</div>


### Ionic Equilibrium

* Total number of cations $n_+$ and anions $n_-$, effective charge of protein $n_p$, then electroneutrality needs $n_+ + n_- = zn_p$.

* Use Donnan approximation which assumes a fixed electrostatic potential of cell $\psi < 0$.


* Cation equilibrium
<div class="tooltip">$k_BT\ln n_+a^3 + e\psi = k_BT\ln n_ba^3 + \Delta\mu_+$
   <span class="tooltiptext">$n_b$ buffer ion concentration, $a$ is microscopic length</span>
</div>

* Anion equilibrium
<div class="tooltip">$k_BT\ln n_-a^3 - e\psi = k_BT\ln n_ba^3 + \Delta\mu_-$


### Osmotic Pressure

* Osmotic pressure inside the cell $\Pi_{in} = k_BT(n_p + n_+ + n_-)$ and that outside $\Pi_{out} = 2n_bk_BT$

* Using Laplace's law one can obtain
<div class="tooltip">$\frac{\Pi_{in}}{k_BT} = \frac{N_p}{V - V_m} + \sqrt{(\frac{zN_p}{V-V_m})^2 + (2n_be^\delta)^2}$
   <span class="tooltiptext">$\delta = \frac{\Delta\mu_+ + \Delta\mu_-}{2k_BT}$</span>
</div>

* When $V \rightarrow V_m$ the osmotic pressure inside the cell reduces to
$\Pi_{in} = \frac{(1 + z)N_pk_BT}{V - V_m}$
	


### Activity and Ion Transport Rates

* Transport rate of ions across $N_{channels}$ of length $h$ is
<div class="tooltip">$r_{channel} = N_{channel}|\delta|/\tau_D$
   <span class="tooltiptext">$\tau_D = \frac{h^2}{D}$</span>
</div>

* Transport rate of ions across $N_{pumps}$ 
<div class="tooltip">$r_{pumps} = rN_{pumps}$
   <span class="tooltiptext">$r$ is pumping rate</span>
</div>

* At steady state $r_{pumps} = r_{channels}$ which gives
<div class="tooltip">$$\delta = -\frac{N_{pumps}}{N_{channels}}r\tau_D$$
   <span class="tooltiptext">Diffusion time $\tau_D$ constant, $N_{pumps}/N_{channels}$ microscopic property</span>
</div>


### Activity and Basal Area

* The response of ion channels depends on membrane+cortical tension which varies with area.

* The number of channels also directly varies with area.

* A simple linear form of this dependence is assumed here
<div class="tooltip">$$\delta = \delta_0 + \delta_1\frac{A - A_s}{A_s} = \delta_0 + \delta_1(\tilde{A} - 1)$$
   <span class="tooltiptext">$A_s = V_s^{2/3}$, where $V$ = volume in suspension, $\tilde{A} = A/A_s$</span>
</div>


### Two Volume Scales

* At the limit $P_{out} \rightarrow \infty$, $V \rightarrow V_m = N_pv$, the "dried" volume.

* In the limit $exp(2\delta) \rightarrow 0$, $V \rightarrow V_0 = V_m + (1 + z)N_p/2n_b$ where $V_0$ is the "desalinated" volume.

* Combining these, define a dimension volume $\tilde{V}$
$$\tilde{V} = \frac{V - V_0}{V_0 - V_m}$$


### Final Result

* Applying hydrostatic equilibrium the following relation is obtained between $\tilde{A}$ and $\tilde{V}$

<div class="tooltip"><div class="testbox">$$e^{2\delta\tilde{A}} = \frac{f(\tilde{V})}{f(\tilde{V_s})}$$</div>
   <span class="tooltiptext">$f(y) = \frac{(1 + z)y^2 + 2zy}{(1 + y)^2}$</span>
</div>


### Model Parameters and Limits

* Measurable volume scales $V_m$, $V_0$, $V_s$.

* Phenomenological parameters $\delta_1$ and $-ze$

* For spread cells $\tilde{V} << 1 \approx \exp(-2|\delta_1|\tilde{A})$

* For $\tilde{A} << 1$, $\tilde{V} - \tilde{V_s} \approx \tilde{A}$



# Comparison with Experiments


## Volume vs Area

![Volume vs. Area](presentation/images/Volume_vs_Area.jpg "Volume vs Area")


## Volume vs Pressure

<img src="presentation/images/Volume_vs_Pressure.jpg" width=60%>



## Conclusions

* The Volume-Area dependence in the model arises due to
  * Activity of ions and channels changing with area.
  * This activity changes the electrochemical potential difference.
  
* Cell shape does not explicitly enter the model due to *constant* number of proteins in the cell. It is captured in the fitted parameter $\delta_1$.


* The assumed linear relation between $\delta_1$ and basal area is valid for *flat cells with comparable apical and basal area.

* The model can be tested by a) Slowly varying buffer pressure and measuring cell volume b) Changing protein content of the cells and measuring cell volume.



## Evaluation and Extension

* Model is simple, analytically tractable and based on experimentally measurable quantities.

* Assumption of fixed nuclear volume may not hold for migrating cells as well as nuclear to cytoplasmic transport.

* No consideration of cell volume regulatory mechanisms dueing cell cycle progression.


* Protein content of cells will change (particularly for aneuploid cells) and due to complex formation which will bring cell shape into play.

* A more fundamental derivation of relationship between activity $\delta_1$ and basal area from microscopic considerations may be necessary.

* For cells interacting in a tissue many of the key assumptions of constant buffer concentration may not hold true.



# Thank You!
