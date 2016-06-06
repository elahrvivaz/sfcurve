/***********************************************************************
* Copyright (c) 2013-2015 Commonwealth Computer Research, Inc.
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Apache License, Version 2.0 which
* accompanies this distribution and is available at
* http://www.opensource.org/licenses/apache2.0.php.
*************************************************************************/
package org.locationtech.sfcurve.zorder

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

object Z3 extends ZObject[Z3Range] {

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

  // Z3 specific methods - we can't generalize these without incurring allocation costs

  /**
    * Create a zrange for this z object
    *
    * @param min min value in range
    * @param max max value in range
    * @return
    */
  def zRange(min: Z3, max: Z3): Z3Range = zRange(min.z, max.z)

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
