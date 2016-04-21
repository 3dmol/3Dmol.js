(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (factory((global.MMTF = global.MMTF || {})));
}(this, function (exports) { 'use strict';

  /**
   * @file msgpack-decode
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * decode binary encoded MessagePack v5 (http://msgpack.org/) data
   * @param  {Uint8Array} buffer - binary encoded MessagePack data
   * @return {Object|Array|String|Number|Boolean|null} decoded Messagepack data
   */
  function decodeMsgpack(buffer) {
    // Loosely based on
    // The MIT License (MIT)
    // Copyright (c) 2013 Tim Caswell <tim@creationix.com>
    // https://github.com/creationix/msgpack-js
    var offset = 0;
    var dataView = new DataView(buffer.buffer);

    /**
     * decode all key-value pairs of a map into an object
     * @param  {Integer} length - number of key-value pairs
     * @return {Object} decoded map
     */
    function map(length) {
      var value = {};
      for (var i = 0; i < length; i++) {
        var key = parse();
        value[key] = parse();
      }
      return value;
    }

    /**
     * decode binary array
     * @param  {Integer} length - number of elements in the array
     * @return {Uint8Array} decoded array
     */
    function bin(length) {
      var value = buffer.subarray(offset, offset + length);
      offset += length;
      return value;
    }

    /**
     * decode string
     * @param  {Integer} length - number string characters
     * @return {String} decoded string
     */
    function str(length) {
      var array = buffer.subarray(offset, offset + length);
      offset += length;
      // limit number of arguments to String.fromCharCode to something
      // browsers can handle, see http://stackoverflow.com/a/22747272
      var chunkSize = 0xffff;
      if(length > chunkSize){
        var c = [];
        for(var i = 0; i < array.length; i += chunkSize) {
          c.push(String.fromCharCode.apply(
            null, array.subarray(i, i + chunkSize)
          ));
        }
        return c.join("");
      }else{
        return String.fromCharCode.apply(null, array);
      }
    }

    /**
     * decode array
     * @param  {Integer} length - number of array elements
     * @return {Array} decoded array
     */
    function array(length) {
      var value = new Array(length);
      for (var i = 0; i < length; i++) {
        value[i] = parse();
      }
      return value;
    }

    /**
     * recursively parse the MessagePack data
     * @return {Object|Array|String|Number|Boolean|null} decoded MessagePack data
     */
    function parse() {
      var type = buffer[offset];
      var value, length, extType;
      // Positive FixInt
      if ((type & 0x80) === 0x00) {
        offset++;
        return type;
      }
      // FixMap
      if ((type & 0xf0) === 0x80) {
        length = type & 0x0f;
        offset++;
        return map(length);
      }
      // FixArray
      if ((type & 0xf0) === 0x90) {
        length = type & 0x0f;
        offset++;
        return array(length);
      }
      // FixStr
      if ((type & 0xe0) === 0xa0) {
        length = type & 0x1f;
        offset++;
        return str(length);
      }
      // Negative FixInt
      if ((type & 0xe0) === 0xe0) {
        value = dataView.getInt8(offset);
        offset++;
        return value;
      }
      switch (type) {
      // nil
      case 0xc0:
        offset++;
        return null;
      // 0xc1: (never used, could be employed for padding)
      // false
      case 0xc2:
        offset++;
        return false;
      // true
      case 0xc3:
        offset++;
        return true;
      // bin 8
      case 0xc4:
        length = dataView.getUint8(offset + 1);
        offset += 2;
        return bin(length);
      // bin 16
      case 0xc5:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return bin(length);
      // bin 32
      case 0xc6:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return bin(length);
      // ext 8
      case 0xc7:
        length = dataView.getUint8(offset + 1);
        extType = dataView.getUint8(offset + 2);
        offset += 3;
        return [extType, bin(length)];
      // ext 16
      case 0xc8:
        length = dataView.getUint16(offset + 1);
        extType = dataView.getUint8(offset + 3);
        offset += 4;
        return [extType, bin(length)];
      // ext 32
      case 0xc9:
        length = dataView.getUint32(offset + 1);
        extType = dataView.getUint8(offset + 5);
        offset += 6;
        return [extType, bin(length)];
      // float 32
      case 0xca:
        value = dataView.getFloat32(offset + 1);
        offset += 5;
        return value;
      // float 64
      case 0xcb:
        value = dataView.getFloat64(offset + 1);
        offset += 9;
        return value;
      // uint8
      case 0xcc:
        value = buffer[offset + 1];
        offset += 2;
        return value;
      // uint 16
      case 0xcd:
        value = dataView.getUint16(offset + 1);
        offset += 3;
        return value;
      // uint 32
      case 0xce:
        value = dataView.getUint32(offset + 1);
        offset += 5;
        return value;
      // uint64
      case 0xcf:
        // FIXME not available/representable in JS
        // largest possible int in JS is 2^53
        // value = dataView.getUint64(offset + 1);
        offset += 9;
        return 0;
      // int 8
      case 0xd0:
        value = dataView.getInt8(offset + 1);
        offset += 2;
        return value;
      // int 16
      case 0xd1:
        value = dataView.getInt16(offset + 1);
        offset += 3;
        return value;
      // int 32
      case 0xd2:
        value = dataView.getInt32(offset + 1);
        offset += 5;
        return value;
      // int 64
      case 0xd3:
        // FIXME not available/representable in JS
        // largest possible int in JS is 2^53
        // value = dataView.getInt64(offset + 1);
        offset += 9;
        return 0;

      // fixext 1
      case 0xd4:
        extType = dataView.getUint8(offset + 1);
        offset += 2;
        return [extType, bin(1)];
      // fixext 2
      case 0xd5:
        extType = dataView.getUint8(offset + 1);
        offset += 2;
        return [extType, bin(2)];
      // fixext 4
      case 0xd6:
        extType = dataView.getUint8(offset + 1);
        offset += 2;
        return [extType, bin(4)];
      // fixext 8
      case 0xd7:
        extType = dataView.getUint8(offset + 1);
        offset += 2;
        return [extType, bin(8)];
      // fixext 16
      case 0xd8:
        extType = dataView.getUint8(offset + 1);
        offset += 2;
        return [extType, bin(16)];
      // str 8
      case 0xd9:
        length = dataView.getUint8(offset + 1);
        offset += 2;
        return str(length);
      // str 16
      case 0xda:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return str(length);
      // str 32
      case 0xdb:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return str(length);
      // array 16
      case 0xdc:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return array(length);
      // array 32
      case 0xdd:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return array(length);
      // map 16:
      case 0xde:
        length = dataView.getUint16(offset + 1);
        offset += 3;
        return map(length);
      // map 32
      case 0xdf:
        length = dataView.getUint32(offset + 1);
        offset += 5;
        return map(length);
      }

      throw new Error("Unknown type 0x" + type.toString(16));
    }

    // start the recursive parsing
    return parse();
  }

  /**
   * @file mmtf-utils
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * get an Uint8Array view on the input array memory
   * @param  {TypedArray} dataArray - input array
   * @return {Uint8Array} new view on the input array memory
   */
  function getUint8View( dataArray ){
      return new Uint8Array(
          dataArray.buffer, dataArray.byteOffset, dataArray.byteLength
      );
  }

  /**
   * get an Int8Array view on the input array memory
   * @param  {TypedArray} dataArray - input array
   * @return {Int8Array} new view on the input array memory
   */
  function getInt8View( dataArray ){
      return new Int8Array(
          dataArray.buffer, dataArray.byteOffset, dataArray.byteLength
      );
  }

  /**
   * get an Int16Array copy of the the input array data
   * @param  {TypedArray} view - input data in big endian format
   * @param  {Int16Array} [dataArray] - pre-allocated output array
   * @return {Int16Array} copy of the input array data
   */
  function getInt16( view, dataArray ){
      var o = view.byteOffset;
      var n = view.byteLength;
      var i, i2, il;
      if( !dataArray ) dataArray = new Int16Array( n / 2 );
      for( i = 0, i2 = 0, il = n / 2; i < il; ++i, i2 += 2 ){
          dataArray[ i ] = view[ i2 ] << 8 ^ view[ i2 + 1 ] << 0;
      }
      return dataArray;
  }

  /**
   * get an Int32Array copy of the the input array data
   * @param  {TypedArray} view - input data in big endian format
   * @param  {Int32Array} [dataArray] - pre-allocated output array
   * @return {Int32Array} copy of the input array data
   */
  function getInt32( view, dataArray ){
      var o = view.byteOffset;
      var n = view.byteLength;
      var i, i4, il;
      if( !dataArray ) dataArray = new Int32Array( n / 4 );
      for( i = 0, i4 = 0, il = n / 4; i < il; ++i, i4 += 4 ){
          dataArray[ i ] = (
              view[ i4     ] << 24 ^ view[ i4 + 1 ] << 16 ^
              view[ i4 + 2 ] <<  8 ^ view[ i4 + 3 ] <<  0
          );
      }
      return dataArray;
  }

  /**
   * get an Int32Array view on the input array memory
   * @param  {TypedArray} dataArray - input array
   * @return {Int32Array} new view on the input array memory
   */
  function getInt32View( dataArray ){
      return new Int32Array(
          dataArray.buffer, dataArray.byteOffset, dataArray.byteLength/4
      );
  }

  /**
   * decode integers into floats using given divisor
   * example:
   *     intArray: [ 12, 34, 543, 687, 2, 0, 4689 ]
   *     divisor: 100
   *     return: [ 0.12, 0.34, 5.43, 6.87, 0.02, 0.00, 46.89 ]
   * @param  {TypedArray|Array} intArray - input array containing integers
   * @param  {Number} divisor - number to devide the integers to obtain floats
   * @param  {Float32Array} [dataArray] - pre-allocated output array
   * @return {Float32Array} decoded array
   */
  function decodeIntegerToFloat( intArray, divisor, dataArray ){
      var n = intArray.length;
      var invDiv = 1/divisor;
      if( !dataArray ) dataArray = new Float32Array( n );
      for( var i = 0; i < n; ++i ){
          // multiply by inverse of the divisor which is faster the division
          dataArray[ i ] = intArray[ i ] * invDiv;
      }
      return dataArray;
  }

  /**
   * perform run-length decoding of input array
   * example:
   *     array: [ 0, 2, 3, 5 ]  // pairs of values and length of a run
   *     return: [ 0, 0, 3, 3, 3, 3, 3 ]
   * @param  {TypedArray|Array} array - run-length encoded input array
   * @param  {TypedArray|Array} [dataArray] - pre-allocated output array
   * @return {TypedArray|Array} decoded array
   */
  function decodeRunLength( array, dataArray ){
      var i, il;
      if( !dataArray ){
          // calculate the length the decoded array will have
          var fullLength = 0;
          for( i = 0, il = array.length; i < il; i+=2 ){
              fullLength += array[ i + 1 ];
          }
          // create a new array of the same type of the input array
          dataArray = new array.constructor( fullLength );
      }
      var dataOffset = 0;
      for( i = 0, il = array.length; i < il; i+=2 ){
          var value = array[ i ];  // value to be repeated
          var length = array[ i + 1 ];  // number of repeats
          for( var j = 0; j < length; ++j ){
              dataArray[ dataOffset ] = value;
              dataOffset += 1;
          }
      }
      return dataArray;
  }

  /**
   * perform delta decoding of the input array
   * by iterativly adding the ith element's value to the i+1th
   * example:
   *     dataArray: [ 0, 2, 1, 2, 1, 1, -4, -2, 9 ]
   *     return: [ 0, 2, 3, 5, 6, 7, 3, 1, 10 ]
   * @param  {TypedArray|Array} dataArray - delta encoded input array
   * @return {TypedArray|Array} decoded array
   */
  function decodeDelta( dataArray ){
      for( var i = 1, il = dataArray.length; i < il; ++i ){
          dataArray[ i ] += dataArray[ i - 1 ];
      }
      return dataArray;
  }

  /**
   * perform split-list delta decoding
   * i.e. the delta values are split between two lists
   * example:
   *     bigArray: [ 200, 3, 100, 2 ]
   *     smallArray: [ 0, 2, -1, -3, 5 ]
   *     return: [ 200, 200, 202, 201, 301, 298, 303 ]
   * @param  {Uint8Array} bigArray - int32 array as bytes in big endian format,
   *                                 pairs of large delta values and number of following
   *                                 small delta values to be read from smallArray
   * @param  {Uint8Array} smallArray - int16 array as bytes in big endian format,
   *                                   small delta values
   * @param  {Int32Array} dataArray - pre-allocated output array
   * @return {Int32Array} decoded array
   */
  function decodeSplitListDelta( bigArray, smallArray, dataArray ){
      var fullLength = ( bigArray.length / 2 ) + smallArray.length;
      if( !dataArray ) dataArray = new Int32Array( fullLength );
      var dataOffset = 0;
      var smallOffset = 0;
      for( var i = 0, il = bigArray.length; i < il; i+=2 ){
          var value = bigArray[ i ];
          var length = bigArray[ i + 1 ];
          dataArray[ dataOffset ] = value;
          if( i !== 0 ){
              dataArray[ dataOffset ] += dataArray[ dataOffset - 1 ];
          }
          dataOffset += 1;
          for( var j = 0; j < length; ++j ){
              dataArray[ dataOffset ] = dataArray[ dataOffset - 1 ] + smallArray[ smallOffset ];
              dataOffset += 1;
              smallOffset += 1;
          }
      }
      return dataArray;
  }

  /**
   * perform split-list delta decoding followed (@see decodeSplitListDelta)
   * by decoding integers into floats using given divisor (@see decodeIntegerToFloat)
   * example:
   *     bigArray: [ 100, 3, -200, 2 ]
   *     smallArray: [ 0, 2, -1, -3, 5 ]
   *     divisor: 100
   *     return: [ 1.00, 1.00, 1.02, 1.01, -0.99, -1.02, -0.97 ]
   * @param  {Uint8Array} bigArray - int32 array as bytes in big endian format,
   *                                 pairs of large delta values and number of following
   *                                 small delta values to be read from smallArray
   * @param  {Uint8Array} smallArray - int16 array as bytes in big endian format,
   *                                   small delta values
   * @param  {Integer} divisor  - number to devide the integers to obtain floats
   * @param  {Float32Array} dataArray - pre-allocated output array
   * @return {Float32Array} decoded array
   */
  function decodeFloatSplitListDelta( bigArray, smallArray, divisor, dataArray ){
      var fullLength = ( bigArray.length / 4 / 2 ) + smallArray.length / 2;
      if( !dataArray ) dataArray = new Float32Array( fullLength );
      var int32View = getInt32View( dataArray );
      var int32 = decodeSplitListDelta(
          getInt32( bigArray ), getInt16( smallArray ), int32View
      );
      return decodeIntegerToFloat( int32, divisor, dataArray );
  }

  /**
   * perform run-length decoding followed (@see decodeRunLength)
   * by decoding integers into floats using given divisor (@see decodeIntegerToFloat)
   * example:
   *     array: [ 320, 3, 100, 2 ]
   *     divisor: 100
   *     return: [ 3.20, 3.20, 3.20, 1.00, 1.00 ]
   * @param  {Uint8Array} array - run-length encoded int32 array as bytes in big endian format
   * @param  {Integer} divisor - number to devide the integers to obtain floats
   * @param  {Float32Array} dataArray - pre-allocated output array
   * @return {Float32Array} decoded array
   */
  function decodeFloatRunLength( bytes, divisor, dataArray ){
      var int32View = dataArray ? getInt32View( dataArray ) : undefined;
      var int32 = decodeRunLength( getInt32( bytes ), int32View );
      return decodeIntegerToFloat( int32, divisor, dataArray );
  }

  /**
   * Decode MMTF fields
   * @param  {Object} inputDict - encoded MMTF data
   * @param  {Object} [params] - decoding parameters
   *  - @param {Array} params.ignoreFields - names of optional fields not to decode
   * @return {Object} mmtfData
   */
  function decodeMmtf( inputDict, params ){

      params = params || {};

      var ignoreFields = params.ignoreFields;

      // helper function to tell if a field should be decoded
      function decodeField( name ){
          return ignoreFields ? ignoreFields.indexOf( name ) === -1 : true;
      }

      // hoisted loop variables
      var i, il;

      // get counts
      var numBonds = inputDict.numBonds || 0;
      var numAtoms = inputDict.numAtoms || 0;
      var numGroups = inputDict.groupTypeList.length / 4;
      var numChains = inputDict.chainIdList.length / 4;
      var numModels = inputDict.chainsPerModel.length;

      // initialize output dict
      var outputDict = {
          numGroups: numGroups,
          numChains: numChains,
          numModels: numModels
      };

      // copy some fields over from the input dict
      [
          "mmtfVersion", "mmtfProducer",
          "unitCell", "spaceGroup", "structureId", "title",
          "depositionDate", "releaseDate",
          "experimentalMethods", "resolution", "rFree", "rWork",
          "bioAssemblyList", "entityList", "groupList",
          "numBonds", "numAtoms",
          "groupsPerChain", "chainsPerModel"
      ].forEach( function( name ){
          if( inputDict[ name ] !== undefined ){
              outputDict[ name ] = inputDict[ name ];
          }
      } );

      //////////////
      // bond data

      // decode inter group bond atom indices, i.e. get int32 array
      var bondAtomListKey = "bondAtomList";
      if( inputDict[ bondAtomListKey ] && decodeField( bondAtomListKey ) ){
          outputDict[ bondAtomListKey ] = getInt32( inputDict[ bondAtomListKey ] );
      }

      // decode inter group bond orders, i.e. get uint8 array
      var bondOrderListKey = "bondOrderList";
      if( inputDict[ bondOrderListKey ] && decodeField( bondOrderListKey ) ){
          outputDict[ bondOrderListKey ] = getUint8View( inputDict[ bondOrderListKey ] );
      }

      //////////////
      // atom data

      // split-list delta & integer decode x, y, z atom coords
      outputDict.xCoordList = decodeFloatSplitListDelta(
          inputDict.xCoordBig, inputDict.xCoordSmall, 1000
      );
      outputDict.yCoordList = decodeFloatSplitListDelta(
          inputDict.yCoordBig, inputDict.yCoordSmall, 1000
      );
      outputDict.zCoordList = decodeFloatSplitListDelta(
          inputDict.zCoordBig, inputDict.zCoordSmall, 1000
      );

      // split-list delta & integer decode b-factors
      var bFactorListKey = "bFactorList";
      var bFactorBigKey = "bFactorBig";
      var bFactorSmallKey = "bFactorSmall";
      if( inputDict[ bFactorBigKey ] && inputDict[ bFactorSmallKey ] && decodeField( bFactorListKey ) ){
          outputDict[ bFactorListKey ] = decodeFloatSplitListDelta(
              inputDict[ bFactorBigKey ], inputDict[ bFactorSmallKey ], 100
          );
      }

      // delta & run-length decode atom ids
      var atomIdListKey = "atomIdList";
      if( inputDict[ atomIdListKey ] && decodeField( atomIdListKey ) ){
          outputDict[ atomIdListKey ] = decodeDelta(
              decodeRunLength( getInt32( inputDict[ atomIdListKey ] ) )
          );
      }

      // run-length decode alternate labels
      var altLocListKey = "altLocList";
      if( inputDict[ altLocListKey ] && decodeField( altLocListKey ) ){
          outputDict[ altLocListKey ] = decodeRunLength(
              getInt32( inputDict[ altLocListKey ] ), new Uint8Array( numAtoms )
          );
      }

      // run-length & integer decode occupancies
      var occupancyListKey = "occupancyList";
      if( inputDict[ occupancyListKey ] && decodeField( occupancyListKey ) ){
          outputDict[ occupancyListKey ] = decodeFloatRunLength( inputDict[ occupancyListKey ], 100 );
      }

      ///////////////
      // group data

      // run-length & delta decode group numbers
      outputDict.groupIdList = decodeDelta(
          decodeRunLength( getInt32( inputDict.groupIdList ) )
      );

      // decode group types, i.e. get int32 array
      outputDict.groupTypeList = getInt32( inputDict.groupTypeList );

      // decode secondary structure, i.e. get int8 view
      var secStructListKey = "secStructList";
      if( inputDict[ secStructListKey ] && decodeField( secStructListKey ) ){
          outputDict[ secStructListKey ] = getInt8View( inputDict[ secStructListKey ] );
      }

      // run-length decode insertion codes
      var insCodeListKey = "insCodeList";
      if( inputDict[ insCodeListKey ] && decodeField( insCodeListKey ) ){
          outputDict[ insCodeListKey ] = decodeRunLength(
              getInt32( inputDict[ insCodeListKey ] ), new Uint8Array( numGroups )
          );
      }

      // run-length & delta decode sequence indices
      var sequenceIndexListKey = "sequenceIndexList";
      if( inputDict[ sequenceIndexListKey ] && decodeField( sequenceIndexListKey ) ){
          outputDict[ sequenceIndexListKey ] = decodeDelta(
              decodeRunLength( getInt32( inputDict[ sequenceIndexListKey ] ) )
          );
      }

      ///////////////
      // chain data

      // decode chain ids, i.e. get int8 view
      outputDict.chainIdList = getUint8View( inputDict.chainIdList );

      // decode chain names, i.e. get int8 view
      var chainNameListKey = "chainNameList";
      if( inputDict[ chainNameListKey ] && decodeField( chainNameListKey ) ){
          outputDict[ chainNameListKey ] = getUint8View( inputDict[ chainNameListKey ] );
      }

      return outputDict;

  }

  /**
   * @file mmtf-iterator
   * @author Alexander Rose <alexander.rose@weirdbyte.de>
   */

  /**
   * Helper methods to loop over MMTF data
   * @class
   * @param {Object} mmtfData - decoded mmtf data object
   */
  function MmtfIterator( mmtfData ){

      var d = mmtfData;

      /**
       * Converts an array of ASCII codes trimming '\0' bytes
       * @param  {Array} charCodeArray - array of ASCII char codes
       * @return {String} '\0' trimmed string
       */
      function fromCharCode( charCodeArray ){
          return String.fromCharCode.apply( null, charCodeArray ).replace(/\0/g, '');
      }

      /**
       * Invokes the callback for each bond
       * @param  {Function} callback(arg0, arg1, arg2) - called for each bond
       *  - @param {Integer} arg0 - first atom index of the bond
       *  - @param {Integer} arg1 - second atom index of the bond
       *  - @param {Integer|null} arg2 - order of the bond
       */
      function eachBond( callback ){
          var i, il;
          // intra group bonds
          var atomOffset = 0;
          for( i = 0, il = d.numGroups; i < il; ++i ){
              var groupData = d.groupList[ d.groupTypeList[ i ] ];
              for( var j = 0, jl = groupData.bondOrderList.length; j < jl; ++j ){
                  callback(
                      atomOffset + groupData.bondAtomList[ j * 2 ],
                      atomOffset + groupData.bondAtomList[ j * 2 + 1 ],
                      groupData.bondOrderList[ j ]
                  );
              }
              atomOffset += groupData.atomNameList.length;
          }
          // inter group bonds
          if( d.bondAtomList ){
              for( i = 0, il = d.bondAtomList.length; i < il; i += 2 ){
                  callback(
                      d.bondAtomList[ i ],
                      d.bondAtomList[ i + 1 ],
                      d.bondOrderList ? d.bondOrderList[ i / 2 ] : null
                  );
              }
          }
      }

      /**
       * Invokes the callback for each atom
       * @param  {Function} callback(arg0, ..., arg9) - called for each atom
       *  - @param {String} arg0 - element
       *  - @param {String} arg1 - atom name
       *  - @param {Integer} arg2 - formal charge
       *  - @param {Float} arg3 - x coordinate
       *  - @param {Float} arg4 - y coordinate
       *  - @param {Float} arg5 - z coordinate
       *  - @param {Float|null} arg6 - b-factor
       *  - @param {Integer|null} arg7 - atom id
       *  - @param {Char|null} arg8 - alternate location label
       *  - @param {Float|null} arg9 - occupancy
       */
      function eachAtom( callback ){
          var atomOffset = 0;
          for( var i = 0, il = d.numGroups; i < il; ++i ){
              var groupData = d.groupList[ d.groupTypeList[ i ] ];
              for( var j = 0, jl = groupData.atomNameList.length; j < jl; ++j ){
                  callback(
                      groupData.elementList[ j ],
                      groupData.atomNameList[ j ],
                      groupData.atomChargeList[ j ],
                      d.xCoordList[ atomOffset ],
                      d.yCoordList[ atomOffset ],
                      d.zCoordList[ atomOffset ],
                      d.bFactorList ? d.bFactorList[ atomOffset ] : null,
                      d.atomIdList ? d.atomIdList[ atomOffset ] : null,
                      d.altLocList ? fromCharCode( [ d.altLocList[ atomOffset ] ] ) : null,
                      d.occupancyList ? d.occupancyList[ atomOffset ] : null
                  );
                  atomOffset += 1;
              }
          }
      }

      /**
       * Invokes the callback for each group
       * @param  {Function} callback(arg0, ..., arg9) - called for each group
       *  - @param {String} arg0 - group name
       *  - @param {Char} arg1 - group single letter code
       *  - @param {String} arg2 - chemical component type
       *  - @param {Integer} arg3 - group id
       *  - @param {Integer} arg4 - group type
       *  - @param {Integer|null} arg5 - secondary structure code
       *  - @param {Char|null} arg6 - insertion code
       *  - @param {Integer|null} arg7 - sequence index
       *  - @param {Integer} arg8 - pointer to data of the group's first atom
       *  - @param {Integer} arg9 - number of atoms in the group
       */
      function eachGroup( callback ){
          var atomOffset = 0;
          for( var i = 0, il = d.numGroups; i < il; ++i ){
              var groupData = d.groupList[ d.groupTypeList[ i ] ];
              var groupAtomCount = groupData.atomNameList.length;
              callback(
                  groupData.groupName,
                  groupData.singleLetterCode,
                  groupData.chemCompType,
                  d.groupIdList[ i ],
                  d.groupTypeList[ i ],
                  d.secStructList ? d.secStructList[ i ] : null,
                  d.insCodeList ? fromCharCode( [ d.insCodeList[ i ] ] ) : null,
                  d.sequenceIndexList ? d.sequenceIndexList[ i ] : null,
                  atomOffset,
                  groupAtomCount
              );
              atomOffset += groupAtomCount;
          }
      }

      /**
       * Invokes the callback for each chain
       * @param  {Function} callback(arg0, ..., arg3) - called for each chain
       *  - @param {String} arg0 - chain id
       *  - @param {String|undefined} arg1 - chain name
       *  - @param {Integer} arg2 - pointer to data of the chain's first group
       *  - @param {Integer} arg3 - number of groups in the chain
       */
      function eachChain( callback ){
          var groupOffset = 0;
          for( var i = 0; i < d.numChains; ++i ){
              var chainGroupCount = d.groupsPerChain[ i ];
              callback(
                  fromCharCode( d.chainIdList.subarray( i, i + 4 ) ),
                  d.chainNameList ? fromCharCode( d.chainNameList.subarray( i, i + 4 ) ) : null,
                  groupOffset,
                  chainGroupCount
              );
              groupOffset += chainGroupCount;
          }
      }

      /**
       * Invokes the callback for each model
       * @param  {Function} callback(arg0, arg1) - called for each model
       *  - @param {Integer} arg0 - pointer to data of the models's first chain
       *  - @param {Integer} arg1 - number of chains in the model
       */
      function eachModel( callback ){
          var chainOffset = 0;
          for( var i = 0; i < d.numModels; ++i ){
              var modelChainCount = d.chainsPerModel[ i ];
              callback(
                  chainOffset,
                  modelChainCount
              );
              chainOffset += modelChainCount;
          }
      }

      // API - bind to instance for public access

      this.eachBond = eachBond;
      this.eachAtom = eachAtom;
      this.eachGroup = eachGroup;
      this.eachChain = eachChain;
      this.eachModel = eachModel;

  }

  /**
   * Decode MMTF fields
   * @param  {Uint8Array|ArrayBuffer|Object} binOrDict - binary MessagePack or encoded MMTF data
   * @param  {Object} [params] - decoding parameters
   *  - @param {Array} params.ignoreFields - names of optional fields not to decode
   * @return {Object} mmtfData
   */
  function decode( binOrDict, params ){
  	// make sure binOrDict is not a plain Arraybuffer
      if( binOrDict instanceof ArrayBuffer ){
          binOrDict = new Uint8Array( binOrDict );
      }

      var inputDict;
      if( binOrDict instanceof Uint8Array ){
          // get dict from msgpack
          inputDict = decodeMsgpack( binOrDict );
      }else{
          // already a dict
          inputDict = binOrDict;
      }

      return decodeMmtf( inputDict, params );
  }

  var Iterator = MmtfIterator;

  exports.decode = decode;
  exports.Iterator = Iterator;

}));