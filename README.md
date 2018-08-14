# FlareNet Particle Penetration Calculator (FPPC)

We developed a Python program "FlareNet Particle Penetration Calculator (FPPC)" to calculate line losses of multi-part sampling lines. The target is to facilitate complicated calculations and help future researchers to expedite these cumbersome computations efficiently. Different mechanism cause particles to be lost en route to the sampling device and now with this program user has an estimation about how much of the specific particle diameter has been lost and may correct their measurement. For calculating particle penetration coefficients these loss mechanisms are considered:

1-	Diffusion<br />
2-	Sedimentation (inclined lines included)<br />
3-	Thermophoresis<br />
4-	Turbulent inertial deposition<br />
5-	Inertial deposition: bend<br />
6-	Inertial deposition: contraction<br />
7-	Sampling probe loss

![alt text](https://raw.githubusercontent.com/keyhanB/FlareNet-Particle-Penetration-Calculator/master/Graph%20Output/PASS3%20Line%20-%20Main%20Graph.jpg)

The line consists of multiple sections, and each section has its properties such as inner diameter, bend angle, gas temperature and flow rate, etc. The user should provide these properties to have the proper calculation. The list of the sections can be reached from the excel input folder. This assumption shows that each part is independent of another and the accuracy of the result depends on the accuracy of the input data.

Please do not hesitate to contact writers if you have any comments or recommendation about the program.
