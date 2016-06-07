/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

/**
 * Represents a rectangle in defined by min and max as two opposing points
 * 
 * @param min: lower-left point 
 * @param max: upper-right point
 */
case class Z2Range(min: Z2, max: Z2) {

  require(min.z <= max.z, s"NOT: ${min.z} <= ${max.z}")

  def mid: Z2 = Z2((max.z + min.z)  >>> 1) // overflow safe mean

  def length: Long = max.z - min.z + 1

  // contains in index space (e.g. the long value)
  def contains(bits: Z2): Boolean = bits.z >= min.z && bits.z <= max.z

  // contains in index space (e.g. the long value)
  def contains(r: Z2Range): Boolean =  contains(r.min) && contains(r.max)

  // overlaps in index space (e.g. the long value)
  def overlaps(r: Z2Range): Boolean = contains(r.min) || contains(r.max)

  // contains in user space - each dimension is contained
  def containsInUserSpace(bits: Z2) = {
    val (x, y) = bits.decode
    x >= min.d0 && x <= max.d0 && y >= min.d1 && y <= max.d1
  }

  // contains in user space - each dimension is contained
  def containsInUserSpace(r: Z2Range): Boolean = containsInUserSpace(r.min) && containsInUserSpace(r.max)

  // overlap in user space - if any dimension overlaps
  def overlapsInUserSpace(r: Z2Range): Boolean =
    overlaps(min.d0, max.d0, r.min.d0, r.max.d0) &&
      overlaps(min.d1, max.d1, r.min.d1, r.max.d1)

  private def overlaps(a1: Int, a2: Int, b1: Int, b2: Int) = math.max(a1, b1) <= math.min(a2, b2)
}
