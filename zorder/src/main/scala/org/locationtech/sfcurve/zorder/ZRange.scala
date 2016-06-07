/*
 * Copyright (c) 2013-2016 Commonwealth Computer Research, Inc.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 */

package org.locationtech.sfcurve.zorder

trait ZRange[ZRangeT <: ZRange[ZRangeT]] {

  require(min <= max, s"Not: $min < $max")

  def min: Long
  def max: Long

  def mid: Long = (min + max) >>> 1 // overflow safe mean

  def length: Long = max - min + 1

  // contains in index space (e.g. the long value)
  def contains(z: Long): Boolean = z >= min && z <= max

  // contains in index space (e.g. the long value)
  def contains(r: ZRangeT): Boolean = contains(r.min) && contains(r.max)

  // overlap in index space (e.g. the long value)
  def overlaps(r: ZRangeT): Boolean = contains(r.min) || contains(r.max)

  // contains in user space - each dimension is contained
  def containsInUserSpace(z: Long): Boolean

  // contains in user space - each dimension is contained
  def containsInUserSpace(r: ZRangeT): Boolean = containsInUserSpace(r.min) && containsInUserSpace(r.max)

  // overlap in user space - if any dimension overlaps
  def overlapsInUserSpace(r: ZRangeT): Boolean
}

object Z2Range {
  def apply(min: Z2, max: Z2)(implicit d: DummyImplicit): Z2Range = Z2Range(min.z, max.z)
}

object Z3Range {
  def apply(min: Z3, max: Z3)(implicit d: DummyImplicit): Z3Range = Z3Range(min.z, max.z)
}

case class Z2Range(min: Long, max: Long) extends ZRange[Z2Range] {

  lazy private val min0 = Z2.decode(min, 0)
  lazy private val min1 = Z2.decode(min, 1)
  lazy private val max0 = Z2.decode(max, 0)
  lazy private val max1 = Z2.decode(max, 1)

  override def containsInUserSpace(z: Long): Boolean = {
    val z0 = Z2.decode(z, 0)
    if (z0 < min0 || z0 > max0) { false } else {
      val z1 = Z2.decode(z, 1)
      z1 >= min1 && z1 <= max1
    }
  }

  override def overlapsInUserSpace(r: Z2Range): Boolean = {
    val rmin0 = Z2.decode(r.min, 0)
    val rmax0 = Z2.decode(r.max, 0)
    if (math.max(min0, rmin0) > math.min(max0, rmax0)) { false } else {
      val rmin1 = Z2.decode(r.min, 1)
      val rmax1 = Z2.decode(r.max, 1)
      math.max(min1, rmin1) <= math.min(max1, rmax1)
    }
  }
}

case class Z3Range(min: Long, max: Long) extends ZRange[Z3Range] {

  lazy private val min0 = Z3.decode(min, 0)
  lazy private val min1 = Z3.decode(min, 1)
  lazy private val min2 = Z3.decode(min, 2)
  lazy private val max0 = Z3.decode(max, 0)
  lazy private val max1 = Z3.decode(max, 1)
  lazy private val max2 = Z3.decode(max, 2)

  override def containsInUserSpace(z: Long): Boolean = {
    val z0 = Z3.decode(z, 0)
    if (z0 < min0 || z0 > max0) { false } else {
      val z1 = Z3.decode(z, 1)
      if (z1 < min1 || z1 > max1) { false } else {
        val z2 = Z3.decode(z, 2)
        z2 >= min2 && z2 <= max2
      }
    }
  }

  override def overlapsInUserSpace(r: Z3Range): Boolean = {
    val rmin0 = Z3.decode(r.min, 0)
    val rmax0 = Z3.decode(r.max, 0)
    if (math.max(min0, rmin0) > math.min(max0, rmax0)) { false } else {
      val rmin1 = Z3.decode(r.min, 1)
      val rmax1 = Z3.decode(r.max, 1)
      if (math.max(min1, rmin1) > math.min(max1, rmax1)) { false } else {
        val rmin2 = Z3.decode(r.min, 2)
        val rmax2 = Z3.decode(r.max, 2)
        math.max(min2, rmin2) <= math.min(max2, rmax2)
      }
    }
  }
}

