/*
    GrandProduct(halo2):

                           (A_i + ß) * (S_i + γ)
        Z_next = Z_prev * –––––––––––––––––––––––
                           (A_i' + ß) * (S_i' + γ)

    wrap_into = 1


    GrandSum(mv):
                              1            m(i)
        Z_next = Z_prev + –––––––––  -  ––––––––––
                           α + f_i       α + t_i


    wrap_into = 0
*/

use ark_ff::{BigInteger, Field, PrimeField};
use ark_std::rand::Rng;
use rand::seq::SliceRandom;

/// Generate small table with values from u64
pub fn small_table<R: Rng, F: Field>(n: usize, rng: &mut R) -> Vec<F> {
    (0..n).map(|_| F::from(rng.gen::<u64>())).collect()
}

pub fn count_number_of_1_bits<F: PrimeField>(x: &[F]) -> u32 {
    let num_of_ones = |x: &F| {
        x.into_bigint()
            .to_bits_be()
            .iter()
            .map(|&b| if b { 1 } else { 0 })
            .sum::<u32>()
    };

    x.iter().map(|xi| num_of_ones(xi)).sum()
}

pub fn mock_grand_polys<F: Field, R: Rng>(
    table: &Vec<F>,
    rng: &mut R,
) -> (Vec<F>, Vec<F>) {
    // prepare
    let mut wtns: Vec<F> = Vec::with_capacity(table.len());
    while wtns.len() < table.len() {
        let random_index = rng.gen_range(0..table.len());
        wtns.push(table[random_index]);
    }

    // shuffled witness 
    let mut wtns_sh = wtns.clone();
    wtns_sh.shuffle(rng);

    // shuffled table
    let mut table_sh = table.clone();
    table_sh.shuffle(rng);

    // multiplicities 
    let mut m: Vec<F> = Vec::with_capacity(table.len());
    while m.len() < table.len() {
        // sample mis to be small numbers 
        let mut mi = rng.gen::<u8>();
        m.push(F::from(mi as u64));
    }

    let beta = F::rand(rng);
    let gamma = F::rand(rng);

    let mut grand_prod: Vec<F> = vec![F::one(); table.len()];
    for i in 0..table.len() - 1 {
        grand_prod[i + 1] = grand_prod[i]
            * (wtns[i] + beta)
            * (table[i] + gamma)
            * (wtns_sh[i] + beta).inverse().unwrap()
            * (table_sh[i] + gamma).inverse().unwrap();
    }

    let mut grand_sum: Vec<F> = vec![F::zero(); table.len()];
    for i in 0..table.len() - 1 {
        grand_sum[i + 1] = grand_sum[i] 
            + (wtns[i] + beta).inverse().unwrap() 
            - m[i] * (table[i] + beta).inverse().unwrap();
    }

    (grand_prod, grand_sum)
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr as F;
    use rand::thread_rng;

    use crate::{small_table, mock_grand_polys, count_number_of_1_bits};

    #[test]
    fn random_sample() {
        let k = 14; 
        let n = 1 << k;
        let mut rng = thread_rng();
        let table = small_table::<_, F>(n, &mut rng);

        let (grand_prod, grand_sum) = mock_grand_polys(&table, &mut rng);

        let ones_prod = count_number_of_1_bits(&grand_prod);
        let ones_sum = count_number_of_1_bits(&grand_sum);

        println!("Number of ones in grand_prod: {}", ones_prod);
        println!("Number of ones in grand_sum: {}", ones_sum);
    }
}
