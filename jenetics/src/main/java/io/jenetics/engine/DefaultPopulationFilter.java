package io.jenetics.engine;

import java.util.function.Predicate;
import java.util.function.Supplier;

import io.jenetics.Gene;
import io.jenetics.Phenotype;
import io.jenetics.util.MSeq;
import io.jenetics.util.Seq;

public class DefaultPopulationFilter<
    G extends Gene<?, G>,
    C extends Comparable<? super C>> implements PopulationFilter<G, C> {
  public DefaultPopulationFilter() {}

  @Override
  public FilterResult<G, C> filter(Predicate<? super Phenotype<G, C>> validator,
                                   Seq<Phenotype<G, C>> population,
                                   Supplier<Phenotype<G, C>> phenotypeSupplier,
                                   long maximalPhenotypeAge,
                                   long generation) {
    int killCount = 0;
    int invalidCount = 0;

    final MSeq<Phenotype<G, C>> pop = MSeq.of(population);
    for (int i = 0, n = pop.size(); i < n; ++i) {
      final Phenotype<G, C> individual = pop.get(i);

      if (!validator.test(individual)) {
        pop.set(i, phenotypeSupplier.get());
        ++invalidCount;
      } else if (individual.getAge(generation) > maximalPhenotypeAge) {
        pop.set(i, phenotypeSupplier.get());
        ++killCount;
      }
    }

    return new FilterResult<>(pop.toISeq(), killCount, invalidCount);
  }
}
