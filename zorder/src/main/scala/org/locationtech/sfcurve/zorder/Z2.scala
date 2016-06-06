/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

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

object Z2 extends ZObject[Z2Range] {

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

  // Z2 specific methods - we can't generalize these without incurring allocation costs

  /**
    * Create a zrange for this z object
    *
    * @param min min value in range
    * @param max max value in range
    * @return
    */
  def zRange(min: Z2, max: Z2): Z2Range = zRange(min.z, max.z)

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
