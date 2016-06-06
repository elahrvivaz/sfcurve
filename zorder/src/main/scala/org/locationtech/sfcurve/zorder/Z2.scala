/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

import org.locationtech.sfcurve.IndexRange
import org.locationtech.sfcurve.zorder.Z3.ZPrefix

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

object Z2 {  
  final val MAX_BITS = 31
  final val MAX_MASK = 0x7fffffff // ignore the sign bit, using it breaks < relationship
  final val MAX_DIM = 2
  final val TOTAL_BITS = MAX_BITS * MAX_DIM

  /** insert 0 between every bit in value. Only first 31 bits can be considered. */
  def split(value: Long): Long = {
    var x: Long = value & MAX_MASK  
    x = (x ^ (x << 32)) & 0x00000000ffffffffL
    x = (x ^ (x << 16)) & 0x0000ffff0000ffffL
    x = (x ^ (x <<  8)) & 0x00ff00ff00ff00ffL // 11111111000000001111111100000000..
    x = (x ^ (x <<  4)) & 0x0f0f0f0f0f0f0f0fL // 1111000011110000
    x = (x ^ (x <<  2)) & 0x3333333333333333L // 11001100..
    x = (x ^ (x <<  1)) & 0x5555555555555555L // 1010...
    x
  }

  /** combine every other bit to form a value. Maximum value is 31 bits. */
  def combine(z: Long): Int = {
    var x = z & 0x5555555555555555L
    x = (x ^ (x >>  1)) & 0x3333333333333333L
    x = (x ^ (x >>  2)) & 0x0f0f0f0f0f0f0f0fL
    x = (x ^ (x >>  4)) & 0x00ff00ff00ff00ffL
    x = (x ^ (x >>  8)) & 0x0000ffff0000ffffL
    x = (x ^ (x >> 16)) & 0x00000000ffffffffL
    x.toInt
  }

  def apply(zvalue: Long): Z2 = new Z2(zvalue)

  /**
   * Bits of x and y will be encoded as ....y1x1y0x0
   */
  def apply(x: Int, y:  Int): Z2 = 
    new Z2(split(x) | split(y) << 1)  

  def unapply(z: Z2): Option[(Int, Int)] = 
    Some(z.decode)
  
  def zdivide(p: Z2, rmin: Z2, rmax: Z2): (Z2, Z2) = {
    val (litmax,bigmin) = zdiv(load, MAX_DIM)(p.z, rmin.z, rmax.z)
    (new Z2(litmax), new Z2(bigmin))
  }
  
  /** Loads either 1000... or 0111... into starting at given bit index of a given dimention */
  def load(target: Long, p: Long, bits: Int, dim: Int): Long = {    
    val mask = ~(Z2.split(MAX_MASK >> (MAX_BITS-bits)) << dim)
    val wiped = target & mask
    wiped | (split(p) << dim)
  }

  /**
    * Calculates the longest common binary prefix between two z longs
    *
    * @return (common prefix, number of bits in common)
    */
  def longestCommonPrefix(values: Long*): ZPrefix = {
    var bitShift = Z2.TOTAL_BITS - Z2.MAX_DIM
    var head = values.head >>> bitShift
    while (values.tail.forall(v => (v >>> bitShift) == head) && bitShift > -1) {
      bitShift -= Z2.MAX_DIM
      head = values.head >>> bitShift
    }
    bitShift += Z2.MAX_DIM // increment back to the last valid value
    ZPrefix(values.head & (Long.MaxValue << bitShift), 64 - bitShift)
  }

  def zranges(zbounds: Array[Z2Range],
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

    def checkValue(prefix: Long, quad: Long): Unit = {
      val min: Long = prefix | (quad << offset) // QR + 000...
      val max: Long = min | (1L << offset) - 1 // QR + 111...
      val octRange = Z2Range(new Z2(min), new Z2(max))

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
    offset -= Z2.MAX_DIM

    var level = 0

    while (level < recurseStop && offset >= 0 && !remaining.isEmpty && ranges.size < rangeStop) {
      val next = remaining.poll
      if (next.eq(levelTerminator)) {
        if (!remaining.isEmpty) {
          level += 1
          offset -= Z2.MAX_DIM
          remaining.add(levelTerminator)
        }
      } else {
        val prefix = next._1
        var quad = 0L
        while (quad < 4) {
          checkValue(prefix, quad)
          quad += 1
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
}
