/*
 * Copyright (c) 2013-2016 Commonwealth Computer Research, Inc.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 */

package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.IndexRange

import scala.collection.mutable.ArrayBuffer

trait ZN[ZRangeT <: ZRange[ZRangeT]] {

  // number of bits used to store each dimension
  def BitsPerDimension: Int

  // number of dimensions
  def Dimensions: Int

  // max value for this z object - can be used to mask another long using &
  def MaxMask: Long

  // has to be lazy to allow for initialization in extending traits
  def TotalBits = BitsPerDimension * Dimensions

  /**
    * Insert (Dimensions - 1) zeros between each bit to create a zvalue from a single dimension.
    * Only the first BitsPerDimension can be considered.
    *
    * @param value value to split
    * @return
    */
  def split(value: Long): Long

  /**
    * Combine every (Dimensions - 1) bits to re-create a single dimension. Opposite of split.
    *
    * @param z value to combine
    * @return
    */
  def combine(z: Long): Int

  /**
    * Create a zrange for this z object
    *
    * @param min min value in range
    * @param max max value in range
    * @return
    */
  private [sfcurve] def zRange(min: Long, max: Long): ZRangeT

  /**
    * Decode a single dimension
    *
    * @param z z value
    * @param dimension dimensino to decode - must be in the range [0, Dimensions)
    * @return
    */
  private [sfcurve] def decode(z: Long, dimension: Int): Int =
    if (dimension == 0) combine(z) else combine(z >> dimension)

  /**
    * Returns (litmax, bigmin) for the given range and point
    *
    * @param p point
    * @param rmin minimum value
    * @param rmax maximum value
    * @return (litmax, bigmin)
    */
  private [sfcurve] def untypedZdivide(p: Long, rmin: Long, rmax: Long): (Long, Long) =
    zdiv(load, Dimensions)(p, rmin, rmax)

  /**
    * Calculates ranges in index space that match any of the input bounds. Uses breadth-first searching to
    * allow a limit on the number of ranges returned.
    *
    * To improve performance, the following decisions have been made:
    *   uses loops instead of foreach/maps
    *   uses java queues instead of scala queues
    *   allocates initial sequences of decent size
    *   sorts once at the end before merging
    *
    * @param zbounds search space
    * @param precision precision to consider, in bits (max 64)
    * @param maxRanges loose cap on the number of ranges to return. A higher number of ranges will have less
    *                  false positives, but require more processing.
    * @param maxRecurse max levels of recursion to apply before stopping
    * @return ranges covering the search space
    */
  def zranges(zbounds: Array[ZRangeT],
              precision: Int = 64,
              maxRanges: Option[Int] = None,
              maxRecurse: Option[Int] = Some(ZN.DefaultRecurse)): Seq[IndexRange]

  /**
   * Cuts Z-Range in two and trims based on user space, can be used to perform augmented binary search
   *
   * @param xd: division point
   * @param inRange: is xd in query range
   */
  private [sfcurve] def untypedCut(r: ZRangeT, xd: Long, inRange: Boolean): List[ZRangeT] = {
    if (r.min == r.max) {
      Nil
    } else if (inRange) {
      if (xd == r.min) { // degenerate case, two nodes min has already been counted
        zRange(r.max, r.max) :: Nil
      } else if (xd == r.max) { // degenerate case, two nodes max has already been counted
        zRange(r.min, r.min) :: Nil
      } else {
        zRange(r.min, xd - 1) :: zRange(xd + 1, r.max) :: Nil
      }
    } else {
      val (litmax, bigmin) = untypedZdivide(xd, r.min, r.max)
      zRange(r.min, litmax) :: zRange(bigmin, r.max) :: Nil
    }
  }

  /**
    * Calculates the longest common binary prefix between two z longs
    *
    * @return (common prefix, number of bits in common)
    */
  def longestCommonPrefix(values: Long*): ZPrefix = {
    var bitShift = TotalBits - Dimensions
    var head = values.head >>> bitShift
    while (values.tail.forall(v => (v >>> bitShift) == head) && bitShift > -1) {
      bitShift -= Dimensions
      head = values.head >>> bitShift
    }
    bitShift += Dimensions // increment back to the last valid value
    ZPrefix(values.head & (Long.MaxValue << bitShift), 64 - bitShift)
  }

  /** Loads either 1000... or 0111... into starting at given bit index of a given dimension */
  private def load(target: Long, p: Long, bits: Int, dim: Int): Long = {
    val mask = ~(split(MaxMask >> (BitsPerDimension - bits)) << dim)
    val wiped = target & mask
    wiped | (split(p) << dim)
  }

  /**
   * Implements the the algorithm defined in: Tropf paper to find:
   * LITMAX: maximum z-index in query range smaller than current point, xd
   * BIGMIN: minimum z-index in query range greater than current point, xd
   *
   * @param load: function that knows how to load bits into appropraite dimension of a z-index
   * @param xd: z-index that is outside of the query range
   * @param rmin: minimum z-index of the query range, inclusive
   * @param rmax: maximum z-index of the query range, inclusive
   * @return (LITMAX, BIGMIN)
   */
  private [zorder] def zdiv(load: (Long, Long, Int, Int) => Long, dims: Int)(xd: Long, rmin: Long, rmax: Long): (Long, Long) = {
    require(rmin < rmax, "min ($rmin) must be less than max $(rmax)")
    var zmin: Long = rmin
    var zmax: Long = rmax
    var bigmin: Long = 0L
    var litmax: Long = 0L

    def bit(x: Long, idx: Int) = {
      ((x & (1L << idx)) >> idx).toInt
    }
    def over(bits: Long)  = 1L << (bits - 1)
    def under(bits: Long) = (1L << (bits - 1)) - 1

    var i = 64
    while (i > 0) {
      i -= 1

      val bits = i/dims+1
      val dim  = i%dims

      ( bit(xd, i), bit(zmin, i), bit(zmax, i) ) match {
        case (0, 0, 0) =>
        // continue

        case (0, 0, 1) =>
          zmax   = load(zmax, under(bits), bits, dim)
          bigmin = load(zmin, over(bits), bits, dim)

        case (0, 1, 0) =>
        // sys.error(s"Not possible, MIN <= MAX, (0, 1, 0)  at index $i")

        case (0, 1, 1) =>
          bigmin = zmin
          return (litmax, bigmin)

        case (1, 0, 0) =>
          litmax = zmax
          return (litmax, bigmin)

        case (1, 0, 1) =>
          litmax = load(zmax, under(bits), bits, dim)
          zmin = load(zmin, over(bits), bits, dim)

        case (1, 1, 0) =>
        // sys.error(s"Not possible, MIN <= MAX, (1, 1, 0) at index $i")

        case (1, 1, 1) =>
        // continue
      }
    }
    (litmax, bigmin)
  }
}

object ZN {
  // indicator that we have searched a full level of the quad/oct tree
  val LevelTerminator = (-1L, -1L)

  val DefaultRecurse = 7
}

case class ZPrefix(prefix: Long, precision: Int) // precision in bits