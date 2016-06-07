/***********************************************************************
* Copyright (c) 2013-2015 Commonwealth Computer Research, Inc.
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Apache License, Version 2.0 which
* accompanies this distribution and is available at
* http://www.opensource.org/licenses/apache2.0.php.
*************************************************************************/
package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.IndexRange

import scala.collection.mutable.ArrayBuffer

class Z3(val z: Long) extends AnyVal {
  import Z3._

  def < (other: Z3) = z < other.z
  def > (other: Z3) = z > other.z
  def >= (other: Z3) = z >= other.z
  def <= (other: Z3) = z <= other.z
  def + (offset: Long) = new Z3(z + offset)
  def - (offset: Long) = new Z3(z - offset)
  def == (other: Z3) = other.z == z

  def d0 = combine(z)
  def d1 = combine(z >> 1)
  def d2 = combine(z >> 2)

  def decode: (Int, Int, Int) = (d0, d1, d2)

  def dim(i: Int): Int = if (i == 0) d0 else if (i == 1) d1 else if (i == 2) d2 else {
    throw new IllegalArgumentException(s"Invalid dimension $i - valid dimensions are 0,1,2")
  }

  def inRange(rmin: Z3, rmax: Z3): Boolean = {
    val (x, y, z) = decode
    x >= rmin.d0 &&
      x <= rmax.d0 &&
      y >= rmin.d1 &&
      y <= rmax.d1 &&
      z >= rmin.d2 &&
      z <= rmax.d2
  }

  def mid(p: Z3): Z3 = {
    val (x, y, z) = decode
    val (px, py, pz) = p.decode
    Z3((x + px) >>> 1, (y + py) >>> 1, (z + pz) >>> 1) // overflow safe mean
  }

  def bitsToString = f"(${z.toBinaryString.toLong}%016d)(${d0.toBinaryString.toLong}%08d," +
      f"${d1.toBinaryString.toLong}%08d,${d2.toBinaryString.toLong}%08d)"

  override def toString = f"$z $decode"
}

object Z3 extends ZN[Z3Range] {

  override val Dimensions = 3
  override val BitsPerDimension = 21
  override val MaxMask = 0x1fffffL

  def apply(zvalue: Long) = new Z3(zvalue)

  /**
   * So this represents the order of the tuple, but the bits will be encoded in reverse order:
   *   ....z1y1x1z0y0x0
   * This is a little confusing.
   */
  def apply(x: Int, y:  Int, z: Int): Z3 = new Z3(split(x) | split(y) << 1 | split(z) << 2)

  def unapply(z: Z3): Option[(Int, Int, Int)] = Some(z.decode)

  override def split(value: Long): Long = {
    var x = value & MaxMask
    x = (x | x << 32) & 0x1f00000000ffffL
    x = (x | x << 16) & 0x1f0000ff0000ffL
    x = (x | x << 8)  & 0x100f00f00f00f00fL
    x = (x | x << 4)  & 0x10c30c30c30c30c3L
    (x | x << 2)      & 0x1249249249249249L
  }

  override def combine(z: Long): Int = {
    var x = z & 0x1249249249249249L
    x = (x ^ (x >>  2)) & 0x10c30c30c30c30c3L
    x = (x ^ (x >>  4)) & 0x100f00f00f00f00fL
    x = (x ^ (x >>  8)) & 0x1f0000ff0000ffL
    x = (x ^ (x >> 16)) & 0x1f00000000ffffL
    x = (x ^ (x >> 32)) & MaxMask
    x.toInt
  }

  override def zRange(min: Long, max: Long): Z3Range = Z3Range(min, max)

  override  def zranges(zbounds: Array[Z3Range],
                        precision: Int = 64,
                        maxRanges: Option[Int] = None,
                        maxRecurse: Option[Int] = Some(ZN.DefaultRecurse)): Seq[IndexRange] = {

    import ZN.LevelTerminator

    // stores our results - initial size of 100 in general saves us some re-allocation
    val ranges = new java.util.ArrayList[IndexRange](100)

    // values remaining to process - initial size of 100 in general saves us some re-allocation
    val remaining = new java.util.ArrayDeque[(Long, Long)](100)

    // calculate the common prefix in the z-values - we start processing with the first diff
    val ZPrefix(commonPrefix, commonBits) = longestCommonPrefix(zbounds.flatMap(b => Seq(b.min, b.max)): _*)

    var offset = 64 - commonBits

    // checks if a range is contained in the search space
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

    // checks if a range overlaps the search space
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

    // checks a single value and either:
    //   eliminates it as out of bounds
    //   adds it to our results as fully matching, or
    //   queues up it's children for further processing
    def checkValue(prefix: Long, quadrant: Long): Unit = {
      val min: Long = prefix | (quadrant << offset) // QR + 000...
      val max: Long = min | (1L << offset) - 1 // QR + 111...
      val quadrantRange = Z3Range(min, max)

      if (isContained(quadrantRange) || offset < 64 - precision) {
        // whole range matches, happy day
        ranges.add(IndexRange(quadrantRange.min, quadrantRange.max, contained = true))
      } else if (overlaps(quadrantRange)) {
        // some portion of this range is excluded
        // queue up each sub-range for processing
        remaining.add((min, max))
      }
    }

    // initial level - we just check the single quadrant
    checkValue(commonPrefix, 0)
    remaining.add(LevelTerminator)
    offset -= Dimensions

    // level of recursion
    var level = 0

    val rangeStop = maxRanges.getOrElse(Int.MaxValue)
    val recurseStop = maxRecurse.getOrElse(ZN.DefaultRecurse)

    while (level < recurseStop && offset >= 0 && !remaining.isEmpty && ranges.size < rangeStop) {
      val next = remaining.poll
      if (next.eq(LevelTerminator)) {
        // we've fully processed a level, increment our state
        if (!remaining.isEmpty) {
          level += 1
          offset -= Dimensions
          remaining.add(LevelTerminator)
        }
      } else {
        val prefix = next._1
        checkValue(prefix, 0L)
        checkValue(prefix, 1L)
        checkValue(prefix, 2L)
        checkValue(prefix, 3L)
        checkValue(prefix, 4L)
        checkValue(prefix, 5L)
        checkValue(prefix, 6L)
        checkValue(prefix, 7L)
      }
    }

    // bottom out and get all the ranges that partially overlapped but we didn't fully process
    while (!remaining.isEmpty) {
      val minMax = remaining.poll
      if (!minMax.eq(LevelTerminator)) {
        ranges.add(IndexRange(minMax._1, minMax._2, contained = false))
      }
    }

    // we've got all our ranges - now reduce them down by merging overlapping values
    ranges.sort(IndexRange.IndexRangeIsOrdered)

    var current = ranges.get(0) // note: should always be at least one range
    val result = ArrayBuffer.empty[IndexRange]
    var i = 1
    while (i < ranges.size()) {
      val range = ranges.get(i)
      if (range.lower <= current.upper + 1) {
        // merge the two ranges
        current = IndexRange(current.lower, math.max(current.upper, range.upper), current.contained && range.contained)
      } else {
        // append the last range and set the current range for future merging
        result.append(current)
        current = range
      }
      i += 1
    }
    // append the last range - there will always be one left that wasn't added
    result.append(current)

    result
  }

  // Z3 specific methods - we can't generalize these without incurring allocation costs

  /**
    * Cuts Z-Range in two and trims based on user space, can be used to perform augmented binary search
    *
    * @param xd: division point
    * @param inRange: is xd in query range
    */
  def cut(r: Z3Range, xd: Z3, inRange: Boolean): List[Z3Range] = untypedCut(r, xd.z, inRange)

  /**
    * Returns (litmax, bigmin) for the given range and point
    *
    * @param p point
    * @param rmin minimum value
    * @param rmax maximum value
    * @return (litmax, bigmin)
    */
  def zdivide(p: Z3, rmin: Z3, rmax: Z3): (Z3, Z3) = {
    val (litmax, bigmin) = untypedZdivide(p.z, rmin.z, rmax.z)
    (Z3(litmax), Z3(bigmin))
  }
}
