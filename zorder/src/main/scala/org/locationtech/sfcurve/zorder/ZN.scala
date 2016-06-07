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


object ZN {

  final val MAX_DIM = 3
  final val MAX_BITS = 21
  final val MAX_MASK = 0x1fffffL
  final val TOTAL_BITS = MAX_BITS * MAX_DIM

  def apply(zvalue: Long) = new Z3(zvalue)

  /** insert 00 between every bit in value. Only first 21 bits can be considered. */
  def split(value: Long): Long = {
    var x = value & MAX_MASK
    x = (x | x << 32) & 0x1f00000000ffffL
    x = (x | x << 16) & 0x1f0000ff0000ffL
    x = (x | x << 8)  & 0x100f00f00f00f00fL
    x = (x | x << 4)  & 0x10c30c30c30c30c3L
    (x | x << 2)      & 0x1249249249249249L
  }

  /** combine every third bit to form a value. Maximum value is 21 bits. */
  def combine(z: Long): Int = {
    var x = z & 0x1249249249249249L
    x = (x ^ (x >>  2)) & 0x10c30c30c30c30c3L
    x = (x ^ (x >>  4)) & 0x100f00f00f00f00fL
    x = (x ^ (x >>  8)) & 0x1f0000ff0000ffL
    x = (x ^ (x >> 16)) & 0x1f00000000ffffL
    x = (x ^ (x >> 32)) & MAX_MASK
    x.toInt
  }

  /**
   * So this represents the order of the tuple, but the bits will be encoded in reverse order:
   *   ....z1y1x1z0y0x0
   * This is a little confusing.
   */
  def apply(x: Int, y:  Int, z: Int): Z3 = {
    new Z3(split(x) | split(y) << 1 | split(z) << 2)
  }

  def unapply(z: Z3): Option[(Int, Int, Int)] = Some(z.decode)

  /**
   * Returns (litmax, bigmin) for the given range and point
   */
  def zdivide(p: Z3, rmin: Z3, rmax: Z3): (Z3, Z3) = {
    val (litmax, bigmin) = zdiv(load, MAX_DIM)(p.z, rmin.z, rmax.z)
    (new Z3(litmax), new Z3(bigmin))
  }

  /** Loads either 1000... or 0111... into starting at given bit index of a given dimension */
  private def load(target: Long, p: Long, bits: Int, dim: Int): Long = {
    val mask = ~(Z3.split(MAX_MASK >> (MAX_BITS-bits)) << dim)
    val wiped = target & mask
    wiped | (split(p) << dim)
  }

  /**
    * Calculates the longest common binary prefix between two z longs
    *
    * @return (common prefix, number of bits in common)
    */
  def longestCommonPrefix(values: Long*): ZPrefix = {
    var bitShift = Z3.TOTAL_BITS - Z3.MAX_DIM
    var head = values.head >>> bitShift
    while (values.tail.forall(v => (v >>> bitShift) == head) && bitShift > -1) {
      bitShift -= Z3.MAX_DIM
      head = values.head >>> bitShift
    }
    bitShift += Z3.MAX_DIM // increment back to the last valid value
    ZPrefix(values.head & (Long.MaxValue << bitShift), 64 - bitShift)
  }

  def zranges(zbounds: Array[Z3Range],
              precision: Int = 64,
              maxRanges: Option[Int] = None,
              maxRecurse: Option[Int] = None): Seq[IndexRange] = {

    // calculate the common prefix in the z-values - we start processing with the first diff
    val ZPrefix(commonPrefix, commonBits) = longestCommonPrefix(zbounds.flatMap(b => Seq(b.min.z, b.max.z)): _*)

    val rangeStop = maxRanges.getOrElse(Int.MaxValue)
    val recurseStop = maxRecurse.getOrElse(7)

    // stores our results
    val ranges = new java.util.ArrayList[IndexRange](100)

    var offset = 64 - commonBits

    val levelTerminator = (-1L, -1L)
    val remaining = new java.util.ArrayDeque[(Long, Long)](100)

    def isContained(range: Z3Range): Boolean = {
      var i = 0
      while (i < zbounds.length) {
        if (zbounds(i).containsInUserSpace(range)) {
          return true
        }
        i += 1
      }
      false
    }

    def overlaps(range: Z3Range): Boolean = {
      var i = 0
      while (i < zbounds.length) {
        if (zbounds(i).overlapsInUserSpace(range)) {
          return true
        }
        i += 1
      }
      false
    }

    def checkValue(prefix: Long, oct: Long): Unit = {
      val min: Long = prefix | (oct << offset) // QR + 000...
      val max: Long = min | (1L << offset) - 1 // QR + 111...
      val octRange = Z3Range(new Z3(min), new Z3(max))

      if (isContained(octRange) || offset < 64 - precision) {
        // whole range matches, happy day
        ranges.add(IndexRange(octRange.min.z, octRange.max.z, contained = true))
      } else if (overlaps(octRange)) {
        // some portion of this range is excluded
        // queue up each sub-range for processing
        remaining.add((min, max))
      }
    }

    // initial level
    checkValue(commonPrefix, 0)
    remaining.add(levelTerminator)
    offset -= Z3.MAX_DIM

    var level = 0

    while (level < recurseStop && offset >= 0 && !remaining.isEmpty && ranges.size < rangeStop) {
      val next = remaining.poll
      if (next.eq(levelTerminator)) {
        if (!remaining.isEmpty) {
          level += 1
          offset -= Z3.MAX_DIM
          remaining.add(levelTerminator)
        }
      } else {
        val prefix = next._1
        var oct = 0L
        while (oct < 8) {
          checkValue(prefix, oct)
          oct += 1
        }
      }
    }

    // bottom out and get all the ranges that partially overlapped but we didn't fully process
    while (!remaining.isEmpty) {
      val (min, max) = remaining.poll
      if (min != -1) {
        ranges.add(IndexRange(min, max, contained = false))
      }
    }

    ranges.sort(IndexRange.IndexRangeIsOrdered)

    var current = ranges.get(0)
    val result = ArrayBuffer.empty[IndexRange]
    var i = 1
    while (i < ranges.size()) {
      val range = ranges.get(i)
      if (range.lower <= current.upper + 1) {
        current = IndexRange(current.lower, math.max(current.upper, range.upper), current.contained && range.contained)
      } else {
        result.append(current)
        current = range
      }
      i += 1
    }
    result.append(current)

    result
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
  def zdiv(load: (Long, Long, Int, Int) => Long, dims: Int)(xd: Long, rmin: Long, rmax: Long): (Long, Long) = {
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

  /**
   * Cuts Z-Range in two and trims based on user space, can be used to perform augmented binary search
   *
   * @param xd: division point
   * @param inRange: is xd in query range
   */
  def cut(r: Z3Range, xd: Z3, inRange: Boolean): List[Z3Range] = {
    if (r.min.z == r.max.z) {
      Nil
    } else if (inRange) {
      if (xd.z == r.min.z)      // degenerate case, two nodes min has already been counted
        Z3Range(r.max, r.max) :: Nil
      else if (xd.z == r.max.z) // degenerate case, two nodes max has already been counted
        Z3Range(r.min, r.min) :: Nil
      else
        Z3Range(r.min, xd - 1) :: Z3Range(xd + 1, r.max) :: Nil
    } else {
      val (litmax, bigmin) = Z3.zdivide(xd, r.min, r.max)
      Z3Range(r.min, litmax) :: Z3Range(bigmin, r.max) :: Nil
    }
  }

  case class ZPrefix(prefix: Long, precision: Int) // precision in bits
}
