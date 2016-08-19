# camp-paleon-project-2016

### PDF Legend

1.  Red = Inverse Distance Weight (1/d)
2.  Green = Inverse Distance Squared Weight (1/d^2)
3.  Blue = Simple Average (1)

### Dispersal Kernals Computed:

1.  25 km
2.  50km
3.  75km
4. 100km
5. 150km
6.  200km

## Calibration model

$$ 
y_i \sim \text{binomial}(n_i, p_i)
$$
$$
p_i = \text{inverse-logit}(r_i)
$$
$$
r_i \sim normal(\beta_0 + \beta_1 \cdot r^*_i, \sigma)
$$
$$
\beta_k \sim \text{normal}(0, 100), k \in [0,1]
$$
$$
\sigma \sim \text{uniform}(0,10)
$$
