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

  private val minz = Z2(min)
  private val maxz = Z2(max)

  override def containsInUserSpace(z: Long): Boolean = {
    val (x, y) = Z2(z).decode
    x >= minz.d0 && x <= maxz.d0 && y >= minz.d1 && y <= maxz.d1
  }

  override def overlapsInUserSpace(r: Z2Range): Boolean = {
    overlaps(minz.d0, maxz.d0, Z2(r.min).d0, Z2(r.max).d0) &&
        overlaps(minz.d1, maxz.d1, Z2(r.min).d1, Z2(r.max).d1)
  }

  private def overlaps(a1: Int, a2: Int, b1: Int, b2: Int): Boolean =
    math.max(a1, b1) <= math.min(a2, b2)
}

case class Z3Range(min: Long, max: Long) extends ZRange[Z3Range] {

  private val minz = Z3(min)
  private val maxz = Z3(max)

  override def containsInUserSpace(zv: Long): Boolean = {
    val (x, y, z) = Z3(zv).decode
    x >= minz.d0 && x <= maxz.d0 && y >= minz.d1 && y <= maxz.d1 && z >= minz.d2 && z <= maxz.d2
  }

  override def overlapsInUserSpace(r: Z3Range): Boolean =
    overlaps(minz.d0, maxz.d0, Z3(r.min).d0, Z3(r.max).d0) &&
      overlaps(minz.d1, maxz.d1, Z3(r.min).d1, Z3(r.max).d1) &&
      overlaps(minz.d2, maxz.d2, Z3(r.min).d2, Z3(r.max).d2)

  private def overlaps(a1: Int, a2: Int, b1: Int, b2: Int) = math.max(a1, b1) <= math.min(a2, b2)
}

