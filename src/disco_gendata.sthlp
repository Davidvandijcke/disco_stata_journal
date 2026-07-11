{smcl}
{* *! version 1.2.1 11jul2026}{...}
{viewerjumpto "Syntax" "disco_gendata##syntax"}{...}
{viewerjumpto "Description" "disco_gendata##description"}{...}
{viewerjumpto "Example" "disco_gendata##example"}{...}
{title:Title}

{phang}
{bf:disco_gendata} {hline 2} Generate a synthetic example dataset for {help disco}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:disco_gendata}


{marker description}{...}
{title:Description}

{pstd}
{cmd:disco_gendata} clears the data in memory and creates an artificial panel for
demonstrating and testing {help disco}: 20 units observed over 20 time periods with
50 individual observations per unit-period. Each unit-period draws outcomes {cmd:y}
from a normal distribution with its own mean and standard deviation, and unit 1
receives a treatment effect that increases with {cmd:y} from period 10 onward. The
resulting variables are {cmd:y} (outcome), {cmd:id} (unit identifier), and
{cmd:time} (time period). The random-number seed is set internally so the dataset
is identical across runs.


{marker example}{...}
{title:Example}

{phang2}{cmd:. disco_gendata}{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(10) ci boots(200)}{p_end}
{phang2}{cmd:. disco_estat summary}{p_end}


{title:Author}

{pstd}
David Van Dijcke, University of Virginia. See {help disco} for the accompanying
article and full package documentation.
