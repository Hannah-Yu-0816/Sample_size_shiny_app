# Sample Size Simulation App

Welcome to the documentation page for the Sample Size Simulation Shiny App.

This tool helps researchers estimate the required sample size for detecting treatment effect heterogeneity 
under different simulation settings.

---

## About the App

- **Purpose**: Support research design by simulating scenarios where treatment effects vary across individuals.
- **Target Users**: Clinical trial designers, biostatisticians, epidemiologists.
- **Model Used**: Causal Forest with Honesty Method.

---

## Background and Motivation

In clinical trials, traditional sample size calculations typically target the overall treatment effect 
(Average Treatment Effect, ATE), assuming that all patients respond similarly to treatment. 

However, in reality, treatment effects can vary substantially across individuals. 
This heterogeneity in treatment response necessitates a shift from focusing solely on the overall effect 
to studying heterogeneous treatment effects (HTE).

Our project specifically addresses the following questions:
- What sample size is required to reliably detect heterogeneous treatment effects (HTE)?
- How can simulation-based methods assist in planning studies that account for treatment effect variability?

Through simulation-based approaches, we aim to provide practical guidance for researchers designing studies 
that seek to uncover treatment effect heterogeneity.

---

## How to Use the App Locally

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/sample-size-simulation.git
   ```
2. Open the project folder in RStudio.
3. Install the required R packages.
4. Launch the app

---

## Online Demo
👉 [Try the app online](https://hannahyu354.shinyapps.io/shinyapp/)

---

## Key Simulation Settings Explained
| **Setting**                        | **Description** | **Instruction**
|:------------------------------------|:----------------|:-----------------|
| Number of covariates (d)           | Total covariates generated | |
| Number of effect modifiers (k)     | Covariates that contribute to the treatment effect | |
| Number of other risk factors (d-k) | Covariates that act as distractors in identifying variables relevant to the treatment effect | |
| $\beta_0$ | | |
| $\gamma_0$ | | |
| $\beta_i$ | | |
| $\gamma$ | | |
| X probabilities                     | Control feature distributions for covariates | |
| sample fraction | proportions of samples used to grow each tree | |
| honesty fraction | proportions of samples used for splitting | |

---

## Output Interpretation

Two indices are defined to show the sample size needed to reveal the heterogeneity: 
heterogeneity captured by effect modifiers and success ratio of pure effect modifier leaf. 

### Heterogeneity Captured by Effect Modifiers

The percentage of heterogeneity captured is calculated as follows:

$$
\frac{\frac{1}{n_1}\sum\limits_{i \in S_1}CATE\left(\prod_{j=1}^{k} X_{ij}=1;\underset{\sim}{X_{-ij}}\right)
\-
\frac{1}{n_0}\sum\limits_{i \in S_0}CATE\left(\prod_{j=1}^{k} X_{ij}=0;\underset{\sim}{X_{-ij}}\right)}
{ground \ truth \ of \ CATE}
\times
100\%
$$

Where:
- $S_1$ is the subset where all selected effect modifiers $X_{ij} = 1$
- $S_0$ is the subset where all selected effect modifiers $X_{ij} = 0$
- The ground truth of CATE is obtained from the data simulation process.

### Success Ratio of Pure Effect Modifier Leaf

The second metric focuses on tree structure: We measure the success rate of trees that split on the effect modifier among the forest.

For this, we look at each individual tree in the forest and classify it into one of three categories: success, partial success, and unsuccess.

- Success: The tree successfully isolates the effect modifier/phenotype group.
- Partial Success: The tree isolates only a subset of the effect modifier/phenotype group.
- Failure: The tree splits on noise and cannot isolate the true subgroup at all.

We then calculate the proportion of success, partial success, and failure across all trees in the forest.

Here is an example with just two variables: $X_1$ is the effect modifier, and $X_2$ is other risk factor. 
The eight trees in the following image represent all possible splitting structures for two binary variables.

<img width="1149" alt="Screenshot 2025-04-28 at 11 48 01 PM" src="https://github.com/user-attachments/assets/0c0eee6a-470b-4450-a348-e09be1fc0552" />

- Trees 1 through 4, and Tree 8, are considered success — they isolate effect modifier $X_1$ cleanly
- Trees 6 and 7 are partial success — they isolate a subset of $X_1$ 
- Tree 5 is a failure — it splits only on noise

For more than one effect modifiers, we define success when the tree can find a leaf that seperate all "on" effect modifiers (i.e., equal to 1). 
This mimics a phenotype-driven mechanism, where the treatment only works if a certain biological profile is fully expressed.

---

## Simulation Process

Here’s how the simulation is structured.

We start by evaluating a grid of sample sizes.
At each sample size, we repeat the simulation 500 times (or the time you choose to run) to ensure results are stable and not due to randomness.

For each simulation run, we follow four steps:
1. Generate data with effect modifiers we got from the user
2. Build a causal forest using that data
3. Compute the two heterogeneity indices we defined earlier — heterogeneity captured and tree success rate
4. Finally, we record the result, and after 500 runs, we take the median as the final estimate for that sample size

To make this efficient, everything is implemented in parallel.

---

## Instructions for Users

### Simulation Controls Panel

The Simulation Controls Panel allows users to define the core structure and settings of the simulation. 
These parameters control how the synthetic data is generated and how the causal forest is constructed and applied. 
Below is a step-by-step explanation of each field:

1. Structure Settings

- Number of Effect Modifiers: This setting controls how many covariates act as treatment effect modifiers.
  - The first k covariates (e.g., $X_1$, $X_2$, ..., $X_k$) will be treated as effect modifiers
    For example:
    - If set to 2 → $X_1$ and $X_2$ are effect modifiers
    - If set to 4 → $X_1$, $X_2$, $X_3$ and $X_4$
  · Maximum allowed: 5
    The input will be capped at 5 to limit complexity
  · Adjusting this value will dynamically update:
    · The data generation formula displayed on the right panel
    · The sidebar to show input fields for each effect modifier's parameters

- Number of Other Risk Factors: Defines how many additional covariates affect baseline risk but not the treatment effect.
  · These covariates are labeled $X_(k+1)$, $X_(k+2)$, ..., continuing from the last effect modifier
    For example:
    · If there are 2 effect modifiers → Other Risk Factors start from $X_3$
    · If there are 5 effect modifiers → Other Risk Factors start from $X_6$
  · No upper limit, but increasing this number will add more covariate to the simulation
  · Updating this value will dynamically:
    · Extend the data generation formula shown on the right panel
    · Add new parameter input fields in the sidebar for each added risk factor

After inputing all the parameters and clicking the "Run Simulation" Button, 

Screenshot

### Parameter Tuning Panel

After selecting the overall sample size for your clinical trials. Then you can come to the second panel: Parameter Tuning. 
In this panel, you will have a chance to tune the sample.fraction and honesty.fraction parameters to see their impact on the two indices.

Screenshot
















