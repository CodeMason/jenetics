package io.jenetics.engine;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.Predicate;
import java.util.function.Supplier;

import io.jenetics.Gene;
import io.jenetics.Phenotype;
import io.jenetics.util.MSeq;
import io.jenetics.util.Seq;

public class ForkJoinPopulationFilter<
    G extends Gene<?, G>,
    C extends Comparable<? super C>> implements PopulationFilter<G, C> {
  private final ForkJoinPool forkJoinPool;
  private final int minTaskSize;

  public ForkJoinPopulationFilter(ForkJoinPool forkJoinPool,
                                  int minTaskSize) {
    this.forkJoinPool = forkJoinPool;
    this.minTaskSize = minTaskSize;
  }

  @Override
  public FilterResult<G, C> filter(Predicate<? super Phenotype<G, C>> validator,
                                   Seq<Phenotype<G, C>> population,
                                   Supplier<Phenotype<G, C>> phenotypeSupplier,
                                   long maximalPhenotypeAge,
                                   long generation) {
    MSeq<Phenotype<G, C>> newPopulation = MSeq.of(population);
    FilterParameters<G, C> filterParameters = new FilterParameters<>(
        validator, population, newPopulation, phenotypeSupplier, maximalPhenotypeAge, generation
    );
    int populationSize = population.size();
    CompletableFuture[] futures = new CompletableFuture[populationSize / minTaskSize + 1];
    for (int i = 0, start = 0; start < populationSize; ++i, start += minTaskSize) {
      int stop = Math.min(populationSize, start + minTaskSize);
      futures[i] = CompletableFuture.runAsync(new FilterAction(filterParameters, start, stop));
    }
    CompletableFuture.allOf(futures).join();
    return new FilterResult<>(
        newPopulation.toISeq(),
        (int)filterParameters.killedCount.sum(),
        (int)filterParameters.invalidCount.sum()
    );
  }


  private static class FilterParameters<
      G extends Gene<?, G>,
      C extends Comparable<? super C>> {
    private final Predicate<? super Phenotype<G, C>> validator;
    private final Seq<Phenotype<G, C>> population;
    private final MSeq<Phenotype<G, C>> newPopulation;
    private final Supplier<Phenotype<G, C>> phenotypeSupplier;
    private final LongAdder killedCount = new LongAdder();
    private final LongAdder invalidCount = new LongAdder();
    private final Object oldPopulationLock = new Object();
    private final Object newPopulationLock = new Object();
    private final long maximalPhenotypeAge;
    private final long generation;

    private FilterParameters(Predicate<? super Phenotype<G, C>> validator,
                             Seq<Phenotype<G, C>> population,
                             MSeq<Phenotype<G, C>> newPopulation,
                             Supplier<Phenotype<G, C>> phenotypeSupplier,
                             long maximalPhenotypeAge,
                             long generation) {
      this.validator = validator;
      this.population = population;
      this.newPopulation = newPopulation;
      this.phenotypeSupplier = phenotypeSupplier;
      this.maximalPhenotypeAge = maximalPhenotypeAge;
      this.generation = generation;
    }
  }

  private class FilterAction implements Runnable {
    private final FilterParameters<G, C> filterParameters;
    private final int start;
    private final int stop;

    private FilterAction(FilterParameters<G, C> filterParameters,
                         int start,
                         int stop) {
      this.filterParameters = filterParameters;
      this.start = start;
      this.stop = stop;
    }

    @Override
    public void run() {
      for (int i = start; i < stop; ++i) {
        final Phenotype<G, C> individual = filterParameters.population.get(i);
        if (!isValid(individual)) {
          writeNewPhenotype(i);
          filterParameters.invalidCount.add(1);
        } else if (individual.getAge(filterParameters.generation) > filterParameters.maximalPhenotypeAge) {
          writeNewPhenotype(i);
          filterParameters.killedCount.add(1);
        }
      }
    }

    private boolean isValid(final Phenotype<G, C> individual) {
      return filterParameters.validator.test(individual);
    }

    private void writeNewPhenotype(int i) {
      Phenotype<G, C> newPhenotype = filterParameters.phenotypeSupplier.get();
      filterParameters.newPopulation.set(i, newPhenotype);
    }
  }
}
