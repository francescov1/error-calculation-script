## Error Calculation Matlab Script

This Matlab script can be used to perform uncertainty calculations and produce sensitivity analysis for calculations using up to 5 variables (algorithm used starts to lose accuracy with more than 5 variables)

To use, run the script and enter in a formula in proper Matlab syntax, as prompted. You will then be prompted to enter in each variable value and its respective error. The code will produce a final answer with error and a sensitivity analysis.

<b>Variable labels cannot contain numbers, but can contain multiple letters, both upper and lower case.</b>

Example equation:</br>
OK: `y=4*h/Wr` </br>
Not OK: `y=4*h2/Wr`
