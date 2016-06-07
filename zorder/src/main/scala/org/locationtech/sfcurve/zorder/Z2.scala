/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.IndexRange

import scala.collection.mutable.ArrayBuffer

class Z2(val z: Long) extends AnyVal {
  import Z2._

  def < (other: Z2) = z < other.z
  def > (other: Z2) = z > other.z

  def + (offset: Long) = new Z2(z + offset)
  def - (offset: Long) = new Z2(z - offset)

  def == (other: Z2) = other.z == z

  def decode: (Int, Int) = (combine(z), combine(z>>1))

  def dim(i: Int) = Z2.combine(z >> i)

  def d0 = dim(0)
  def d1 = dim(1)

  def mid(p: Z2): Z2 = {
    val (x, y) = decode
    val (px, py) = p.decode
    Z2((x + px) >>> 1, (y + py) >>> 1) // overflow safe mean
  }

  def bitsToString = f"(${z.toBinaryString}%16s)(${dim(0).toBinaryString}%8s,${dim(1).toBinaryString}%8s)"
  override def toString = f"$z $decode"
}

object Z2 extends ZN[Z2Range] {

  override val Dimensions = 2
  override val BitsPerDimension = 31
  override val MaxMask = 0x7fffffffL // ignore the sign bit, using it breaks < relationship

  def apply(zvalue: Long): Z2 = new Z2(zvalue)

  /**
    * Bits of x and y will be encoded as ....y1x1y0x0
    */
  def apply(x: Int, y:  Int): Z2 = new Z2(split(x) | split(y) << 1)

  def unapply(z: Z2): Option[(Int, Int)] = Some(z.decode)

  override def split(value: Long): Long = {
    var x: Long = value & MaxMask
    x = (x ^ (x << 32)) & 0x00000000ffffffffL
    x = (x ^ (x << 16)) & 0x0000ffff0000ffffL
    x = (x ^ (x <<  8)) & 0x00ff00ff00ff00ffL // 11111111000000001111111100000000..
    x = (x ^ (x <<  4)) & 0x0f0f0f0f0f0f0f0fL // 1111000011110000
    x = (x ^ (x <<  2)) & 0x3333333333333333L // 11001100..
    x = (x ^ (x <<  1)) & 0x5555555555555555L // 1010...
    x
  }

  override def combine(z: Long): Int = {
    var x = z & 0x5555555555555555L
    x = (x ^ (x >>  1)) & 0x3333333333333333L
    x = (x ^ (x >>  2)) & 0x0f0f0f0f0f0f0f0fL
    x = (x ^ (x >>  4)) & 0x00ff00ff00ff00ffL
    x = (x ^ (x >>  8)) & 0x0000ffff0000ffffL
    x = (x ^ (x >> 16)) & 0x00000000ffffffffL
    x.toInt
  }

  override def zRange(min: Long, max: Long): Z2Range = Z2Range(min, max)

  override  def zranges(zbounds: Array[Z2Range],
                        precision: Int = 64,
                        maxRanges: Option[Int],
                        maxRecurse: Option[Int]): Seq[IndexRange] = {

    import ZN.LevelTerminator

    // stores our results - initial size of 100 in general saves us some re-allocation
    val ranges = new java.util.ArrayList[IndexRange](100)

    // values remaining to process - initial size of 100 in general saves us some re-allocation
    val remaining = new java.util.ArrayDeque[(Long, Long)](100)

    // calculate the common prefix in the z-values - we start processing with the first diff
    val ZPrefix(commonPrefix, commonBits) = longestCommonPrefix(zbounds.flatMap(b => Seq(b.min, b.max)): _*)

    var offset = 64 - commonBits

    // checks if a range is contained in the search space
    def isContained(range: Z2Range): Boolean = {
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
    def overlaps(range: Z2Range): Boolean = {
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
      val quadrantRange = Z2Range(min, max)

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

  // Z2 specific methods - we can't generalize these without incurring allocation costs

  /**
    * Create a zrange for this z object
    *
    * @param min min value in range
    * @param max max value in range
    * @return
    */
  def z2Range(min: Z2, max: Z2): Z2Range = Z2Range(min.z, max.z)

  /**
    * Cuts Z-Range in two and trims based on user space, can be used to perform augmented binary search
    *
    * @param xd: division point
    * @param inRange: is xd in query range
    */
  def cut(r: Z2Range, xd: Z2, inRange: Boolean): List[Z2Range] = untypedCut(r, xd.z, inRange)

  /**
    * Returns (litmax, bigmin) for the given range and point
    *
    * @param p point
    * @param rmin minimum value
    * @param rmax maximum value
    * @return (litmax, bigmin)
    */
  def zdivide(p: Z2, rmin: Z2, rmax: Z2): (Z2, Z2) = {
    val (litmax, bigmin) = untypedZdivide(p.z, rmin.z, rmax.z)
    (Z2(litmax), Z2(bigmin))
  }
}
